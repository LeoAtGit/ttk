#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <Propagation.h>

#if(defined(__GNUC__) && !defined(__clang__))
#include <parallel/algorithm>
#endif

// for numerical perturbation
#include <boost/math/special_functions/next.hpp>

typedef ttk::SimplexId ttkInt;

std::string toFixed(const float& number, const int precision = 2){
    std::stringstream vFraction;
    vFraction << std::fixed << std::setprecision(precision) << number;
    return vFraction.str();
};

template<typename idType>
std::string toFixed(const idType& number0, const idType& number1, const int precision = 2){
    return toFixed( ((float)number0)/((float)number1), precision );
};

namespace ttk {

    class LTSimplification : virtual public Debug {

        public:

            LTSimplification(){
                this->setDebugMsgPrefix("LTS"); // inherited from Debug: prefix will be printed at the beginning of every msg
            };
            ~LTSimplification(){};

            int PreconditionTriangulation(
                ttk::Triangulation* triangulation
            ) const {
                triangulation->preconditionVertexNeighbors();
                return 1;
            }

            template<typename T,typename idType>
            int computeGlobalOrder(
                idType* outputOrder,
                std::vector<std::tuple<T,idType,idType>>& sortedIndices,

                const T* rank1,
                const idType* rank2,
                const idType& nVertices,

                std::tuple<idType,idType,idType>* sortedIndicesCopy = nullptr,
                idType* sortedIndicesCopyII = nullptr
            ) const {
                ttk::Timer timer;

                // init tuples
                #pragma omp parallel for num_threads(this->threadNumber_)
                for(idType i=0; i<nVertices; i++){
                    auto& t = sortedIndices[i];
                    std::get<0>(t) = rank1[i];
                    std::get<1>(t) = rank2[i];
                    std::get<2>(t) = i;
                }

                this->printMsg(
                    "Computing global orders",
                    0.2, timer.getElapsedTime(), this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                #ifdef TTK_ENABLE_OPENMP
                    #ifdef __clang__
                        this->printWrn("Caution, outside GCC, sequential sort");
                        std::sort(sortedIndices.begin(), sortedIndices.end());
                    #else
                        omp_set_num_threads(this->threadNumber_);
                        __gnu_parallel::sort(sortedIndices.begin(), sortedIndices.end());
                        omp_set_num_threads(1);
                    #endif
                #else
                    this->printWrn("Caution, outside GCC, sequential sort");
                    std::sort(sortedIndices.begin(), sortedIndices.end());
                #endif

                this->printMsg(
                    "Computing global orders",
                    0.8, timer.getElapsedTime(), this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                if(sortedIndicesCopy){
                    #pragma omp parallel for num_threads(this->threadNumber_)
                    for(idType i=0; i<nVertices; i++){
                        const idType& v = std::get<2>(sortedIndices[i]);
                        outputOrder[v] = i;
                        std::get<2>(sortedIndicesCopy[i]) = v;
                    }
                } else if(sortedIndicesCopyII) {
                    #pragma omp parallel for num_threads(this->threadNumber_)
                    for(idType i=0; i<nVertices; i++){
                        const idType& v = std::get<2>(sortedIndices[i]);
                        outputOrder[v] = i;
                        sortedIndicesCopyII[i] = v;
                    }
                } else {
                    #pragma omp parallel for num_threads(this->threadNumber_)
                    for(idType i=0; i<nVertices; i++)
                        outputOrder[std::get<2>(sortedIndices[i])] = i;
                }

                this->printMsg(
                    "Computing global orders",
                    1,
                    timer.getElapsedTime(),
                    this->threadNumber_
                );

                return 1;
            }

            template<typename idType, typename dataType>
            int enforceAuthorizedExtrema(
                idType* orders,

                const ttk::Triangulation* triangulation,
                const dataType* scalars,
                const idType* authorizedExtremaIndices,
                const idType& nAuthorizedExtremaIndices
            ) const {
                ttk::Timer timer;
                this->printMsg(
                    "Enforcing authorized extrema",
                    0,0,this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                const idType nVertices = triangulation->getNumberOfVertices();

                idType nEnforcedExtrema = 0;

                #pragma omp parallel for reduction(+:nEnforcedExtrema) num_threads(this->threadNumber_)
                for(idType i=0; i<nAuthorizedExtremaIndices; i++){
                    const idType& v = authorizedExtremaIndices[i];

                    const dataType& vScalar = scalars[v];
                    const idType& vOrder = orders[v];

                    // check if v has larger neighbors
                    bool hasLargerNeighborInScalar = false;
                    bool hasLargerNeighborInOrder = false;
                    bool hasSmallerNeighborInScalar = false;
                    bool hasSmallerNeighborInOrder = false;

                    idType nNeighbors = triangulation->getVertexNeighborNumber( v );
                    for(idType n=0; n<nNeighbors; n++){
                        idType u;
                        triangulation->getVertexNeighbor(v,n,u);

                        const dataType& uScalar = scalars[u];
                        const idType& uOrder = orders[u];

                        if( uOrder<vOrder )
                            hasSmallerNeighborInOrder = true;
                        else
                            hasLargerNeighborInOrder = true;

                        if(uScalar<vScalar || (uScalar==vScalar && u<v))
                            hasSmallerNeighborInScalar = true;
                        else
                            hasLargerNeighborInScalar = true;
                    }

                    // if v was a maximum but is not anymore
                    if(!hasLargerNeighborInScalar && hasLargerNeighborInOrder){
                        // force it to be a mximum in the order
                        orders[v] = nVertices;
                        nEnforcedExtrema++;
                    }

                    // if v was a minimum but is not anymore
                    if(!hasSmallerNeighborInScalar && hasSmallerNeighborInOrder) {
                        // force it to be a mximum in the order
                        orders[v] = -nVertices;
                        nEnforcedExtrema++;
                    }
                }

                this->printMsg(
                    "Enforcing authorized extrema ("+std::to_string(nEnforcedExtrema)+")",
                    1,timer.getElapsedTime(),this->threadNumber_
                );

                return 1;
            }

            template<typename dataType, typename idType>
            int computeNumericalPerturbation(
                dataType* outputScalars,

                const idType* orders,
                const std::vector<std::tuple<idType,idType,idType>>& sortedIndices,
                const int sortDirection
            ) const {
                ttk::Timer timer;
                this->printMsg(
                    "Applying numerical perturbation",
                    0,0,this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                const idType nVertices = sortedIndices.size();
                if(sortDirection>0){
                    for(idType i=1; i<nVertices; i++){
                        const idType& v0 = std::get<2>(sortedIndices[i-1]);
                        const idType& v1 = std::get<2>(sortedIndices[i]);
                        if(outputScalars[v0]>=outputScalars[v1])
                            outputScalars[v1] = boost::math::float_next(outputScalars[v0]);
                    }
                } else if(sortDirection<0) {
                    for(idType i=nVertices-1; i>1; i--){
                        const idType& v1 = std::get<2>(sortedIndices[i-1]);
                        const idType& v0 = std::get<2>(sortedIndices[i]);
                        if(outputScalars[v0]>=outputScalars[v1])
                            outputScalars[v1] = boost::math::float_next(outputScalars[v0]);
                    }
                }

                this->printMsg(
                    "Applying numerical perturbation",
                    1,timer.getElapsedTime(),this->threadNumber_
                );

                return 1;
            }

            template<typename idType>
            int flattenOrders(
                idType* outputOrder,

                const std::vector<Propagation<idType>*>& masterPropagations,
                const idType* inputOrders,
                const idType& nVertices
            ) const {
                ttk::Timer timer;
                this->printMsg("Flattening orders",0,0,this->threadNumber_,debug::LineMode::REPLACE);

                const idType nMasterPropagations = masterPropagations.size();

                // use region mask as temporary array
                #pragma omp parallel for num_threads(this->threadNumber_)
                for(idType v=0; v<nVertices; v++)
                    outputOrder[v] = inputOrders[v];

                // flatten regions to order of last encountered saddles
                #pragma omp parallel for num_threads(this->threadNumber_)
                for(idType p=0; p<nMasterPropagations; p++){
                    const auto* propagation = masterPropagations[p];
                    for(const auto& v : propagation->region)
                        outputOrder[v] = inputOrders[propagation->lastEncounteredCriticalPoint];
                }

                this->printMsg("Flattening orders",1,timer.getElapsedTime(),this->threadNumber_);

                return 1;
            }

            template<typename dataType, typename idType>
            int flattenScalars(
                dataType* scalars,

                const std::vector<Propagation<idType>*>& masterPropagations
            ) const {
                ttk::Timer timer;
                this->printMsg("Flattening scalar field",0,0,this->threadNumber_,debug::LineMode::REPLACE);

                const idType nMasterPropagations = masterPropagations.size();

                #pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
                for(idType p=0; p<nMasterPropagations; p++){
                    const auto& propagation = *masterPropagations[p];
                    const idType s = propagation.lastEncounteredCriticalPoint;
                    const dataType sScalar = scalars[s];
                    for(auto v : propagation.region)
                        scalars[v] = sScalar;
                }

                this->printMsg("Flattening scalar field",1,timer.getElapsedTime(),this->threadNumber_);

                return 1;
            }

            template<typename idType>
            int invertField(
                idType* outputOrder,
                idType* inputOrders,

                const idType& nVertices,

                idType* sortedIndices = nullptr
            ) const {
                ttk::Timer timer;
                this->printMsg("Inverting fields",0,0,this->threadNumber_,debug::LineMode::REPLACE);

                if(!sortedIndices){
                    #pragma omp parallel for num_threads(this->threadNumber_)
                    for(idType v=0; v<nVertices; v++){
                        idType& outputOrderV = outputOrder[v];
                        inputOrders[v] = -outputOrderV;
                        outputOrderV = -outputOrderV;
                    }
                } else {
                    #pragma omp parallel for num_threads(this->threadNumber_)
                    for(idType v=0; v<nVertices; v++){
                        idType& outputOrderV = outputOrder[v];
                        inputOrders[v] = -outputOrderV;
                        outputOrderV = -outputOrderV;
                    }
                    #pragma omp parallel for num_threads(this->threadNumber_)
                    for(idType v=0; v<nVertices/2; v++){
                        std::swap(sortedIndices[v], sortedIndices[nVertices-v-1]);
                    }
                }

                this->printMsg("Inverting fields",1,timer.getElapsedTime(),this->threadNumber_);

                return 1;
            }

            template<typename idType>
            int detectUnauthorizedMaxima(
                std::vector<idType>& unauthorizedMaxima,
                idType* authorizationMask,

                const ttk::Triangulation* triangulation,
                const idType* inputOrders,
                const idType* authorizedExtremaIndices,
                const idType& nAuthorizedExtremaIndices
            ) const {

                ttk::Timer timer;
                this->printMsg(
                    "Detecting unauthorized maxima",
                    0,0,this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                const idType nVertices = triangulation->getNumberOfVertices();

                // make room for the maximal number of maxima
                unauthorizedMaxima.resize(nVertices);

                // a synchronized write index used to store discarded maxima
                idType maximaWriteIndex=0;

                // init preservation mask
                #pragma omp parallel for num_threads(this->threadNumber_)
                for(idType i=0; i<nVertices; i++)
                    authorizationMask[i]=-1;

                #pragma omp parallel for num_threads(this->threadNumber_)
                for(idType i=0; i<nAuthorizedExtremaIndices; i++)
                    authorizationMask[authorizedExtremaIndices[i]]=-2;

                // find discareded maxima
                #pragma omp parallel for num_threads(this->threadNumber_)
                for(idType v=0; v<nVertices; v++){

                    // if v needs to be preserved then skip
                    if(authorizationMask[v]==-2)
                        continue;

                    // check if v has larger neighbors
                    bool hasLargerNeighbor = false;
                    const idType& vOrder = inputOrders[v];
                    idType nNeighbors = triangulation->getVertexNeighborNumber( v );
                    for(idType n=0; n<nNeighbors; n++){
                        idType u;
                        triangulation->getVertexNeighbor(v,n,u);
                        if( vOrder<inputOrders[u] ){
                            hasLargerNeighbor = true;
                            break;
                        }
                    }

                    // if v has larger neighbors then v can not be maximum
                    if(hasLargerNeighbor)
                        continue;

                    // get local write index for this thread
                    idType localWriteIndex = 0;
                    #pragma omp atomic capture
                    localWriteIndex = maximaWriteIndex++;

                    // write maximum index
                    unauthorizedMaxima[localWriteIndex] = v;
                }

                // resize to the actual number of discarded maxima
                unauthorizedMaxima.resize(maximaWriteIndex);
                this->printMsg(
                    "Detecting unauthorized maxima ("+std::to_string(maximaWriteIndex)+"|"+std::to_string(nVertices)+")",
                    1,timer.getElapsedTime(),
                    this->threadNumber_
                );

                return 1;
            }

            template<typename idType>
            int detectMaxima(
                std::vector<idType>& maxima,

                const ttk::Triangulation* triangulation,
                const idType* orders
            ) const {

                ttk::Timer timer;
                this->printMsg(
                    "Detecting maxima",
                    0,0,this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                const idType nVertices = triangulation->getNumberOfVertices();

                // make room for the maximal number of maxima
                maxima.resize(nVertices);

                // a synchronized write index used to store discarded maxima
                idType maximaWriteIndex=0;

                // find discareded maxima
                #pragma omp parallel for num_threads(this->threadNumber_)
                for(idType v=0; v<nVertices; v++){
                    // check if v has larger neighbors
                    bool hasLargerNeighbor = false;
                    const idType& vOrder = orders[v];
                    idType nNeighbors = triangulation->getVertexNeighborNumber( v );
                    for(idType n=0; n<nNeighbors; n++){
                        idType u;
                        triangulation->getVertexNeighbor(v,n,u);
                        if( vOrder<orders[u] ){
                            hasLargerNeighbor = true;
                            break;
                        }
                    }

                    // if v has larger neighbors then v can not be maximum
                    if(hasLargerNeighbor)
                        continue;

                    // get local write index for this thread
                    idType localWriteIndex = 0;
                    #pragma omp atomic capture
                    localWriteIndex = maximaWriteIndex++;

                    // write maximum index
                    maxima[localWriteIndex] = v;
                }

                // resize to the actual number of discarded maxima
                maxima.resize(maximaWriteIndex);

                this->printMsg(
                    "Detecting maxima ("+std::to_string(maximaWriteIndex)+"|"+std::to_string(nVertices)+")",
                    1,timer.getElapsedTime(),this->threadNumber_
                );

                return 1;
            }

            template<typename idType>
            int sortMaxima(
                std::vector<idType>& maxima,

                const idType* orders
            ) const {
                ttk::Timer timer;
                this->printMsg(
                    "Sorting maxima ("+std::to_string(maxima.size())+")",
                    0,0,this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                struct Comperator{
                    const idType* orders;
                    Comperator(const idType* newOrders) : orders(newOrders){};
                    inline bool operator ()(const idType& i, const idType& j){
                        return this->orders[i]<this->orders[j];
                    }
                };

                #ifdef TTK_ENABLE_OPENMP
                    #ifdef __clang__
                        this->printWrn("Caution, outside GCC, sequential sort");
                        std::sort(maxima.begin(), maxima.end(), Comperator(orders));
                    #else
                        omp_set_num_threads(this->threadNumber_);
                        __gnu_parallel::sort(maxima.begin(), maxima.end(), Comperator(orders));
                        omp_set_num_threads(1);
                    #endif
                #else
                    this->printWrn("Caution, outside GCC, sequential sort");
                    std::sort(maxima.begin(), maxima.end(), Comperator(orders));
                #endif

                this->printMsg(
                    "Sorting maxima ("+std::to_string(maxima.size())+")",
                    1,timer.getElapsedTime(),this->threadNumber_
                );

                return 1;
            }
            template<typename idType>
            int computeRegion(
                idType* regionMask,
                Propagation<idType>** propagationMask,
                Propagation<idType>* propagation,

                const ttk::Triangulation* triangulation
            ) const {

                // this->printErr("Region size: "+std::to_string(propagation->regionSize)+" "+std::to_string(propagation->extremumIndex)+" -> "+std::to_string(propagation->lastEncounteredCriticalPoint));

                const idType& extremumIndex = propagation->extremumIndex;

                // collect region
                auto& region = propagation->region;

                region.resize(propagation->regionSize);
                idType regionIndex = 0;
                {
                    std::vector<idType> queue(propagation->regionSize,-1);
                    idType queueIndex = 0;
                    {
                        queue[queueIndex++] = extremumIndex;
                        regionMask[extremumIndex] = extremumIndex;

                        if(propagationMask[extremumIndex]->find()!=propagation){
                            this->printErr("Unexpected Configuration");
                        }
                    }

                    while(queueIndex>0){
                        const idType v = queue[--queueIndex];

                        region[regionIndex++] = v;

                        idType nNeighbors = triangulation->getVertexNeighborNumber(v);
                        for(idType n=0; n<nNeighbors; n++){
                            idType u;
                            triangulation->getVertexNeighbor(v,n,u);
                            if(regionMask[u]!=extremumIndex && propagationMask[u]!=nullptr && propagationMask[u]->find()==propagation){
                                queue[queueIndex++]=u;
                                regionMask[u] = extremumIndex;
                            }
                        }
                    }
                }

                if(regionIndex!=propagation->regionSize){
                    this->printErr("Region size incorrect: "+std::to_string(regionIndex)+ " "+std::to_string(propagation->regionSize));
                    return 0;
                }

                return 1;
            }

            template<typename idType>
            int computeRegions(
                idType* regionMask,
                Propagation<idType>** propagationMask,
                std::vector<Propagation<idType>*>& propagations,

                const ttk::Triangulation* triangulation
            ) const {

                const idType nPropagations = propagations.size();
                const idType nVertices = triangulation->getNumberOfVertices();

                ttk::Timer timer;
                this->printMsg(
                    "Computing regions ("+std::to_string(nPropagations)+")",
                    0, 0, this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                int status = 1;

                #pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
                for(idType p=0; p<nPropagations; p++){
                    int localStatus = this->computeRegion<idType>(
                        regionMask,
                        propagationMask,
                        propagations[p],

                        triangulation
                    );
                    if(!localStatus)
                        status = 0;
                }
                if(!status) return 0;


                if(this->debugLevel_<4){
                    this->printMsg(
                        "Computing regions ("+std::to_string(nPropagations)+")",
                        1, timer.getElapsedTime(), this->threadNumber_
                    );
                } else {

                    idType min = propagations[0]->regionSize;
                    idType max = min;
                    idType avg = 0;

                    for(idType p=0; p<nPropagations; p++){
                        const auto propagation = propagations[p];
                        if(min>propagation->regionSize)
                            min=propagation->regionSize;
                        if(max<propagation->regionSize)
                            max=propagation->regionSize;
                        avg += propagation->regionSize;
                    }

                    avg /= nPropagations;

                    this->printMsg(
                        "Computing regions ("+std::to_string(nPropagations)
                            + "|" + toFixed(min,nVertices)
                            + "|" + toFixed(avg,nVertices)
                            + "|" + toFixed(max,nVertices)
                        +")",
                        1, timer.getElapsedTime(), this->threadNumber_
                    );
                }

                return 1;
            }

            template<typename idType>
            inline int getSaddlePropagations(
                std::vector<Propagation<idType>*>& saddlePropagations,
                Propagation<idType>*& masterPropagation,
                Propagation<idType>** propagationMask,

                const ttk::Triangulation* triangulation,
                const idType& saddleIndex,
                const idType& saddleOrder,
                const idType& nNeighbors,
                const idType* orders
            ) const {
                saddlePropagations.resize(nNeighbors);

                for(idType n=0; n<nNeighbors; n++){
                    idType u;
                    triangulation->getVertexNeighbor(saddleIndex,n,u);

                    Propagation<idType>* propagation = nullptr;

                    if(saddleOrder<orders[u]){
                        propagation = propagationMask[u]->find();
                        if(orders[masterPropagation->extremumIndex]<orders[propagation->extremumIndex]){
                            masterPropagation = propagation;
                        }
                    }

                    saddlePropagations[n] = propagation;
                }

                return 1;
            }

            template<typename idType>
            inline int mergeSaddlePropagations(
                std::vector<Propagation<idType>*>& saddlePropagations,
                Propagation<idType>*& masterPropagation
            ) const {
                for(auto* saddlePropagation : saddlePropagations){
                    if(saddlePropagation!=nullptr && saddlePropagation->find()!=masterPropagation){
                        Propagation<idType>::unify2(
                            masterPropagation,
                            saddlePropagation
                        );
                    }
                }

                return 1;
            }

            template<typename idType>
            int computePropagation(
                idType* saddleMask, // used here to store registered larger vertices
                idType* queueMask, // used to mark vertices that have already been added to the queue by this thread
                Propagation<idType>** propagationMask,
                Propagation<idType>& propagation,

                const ttk::Triangulation* triangulation,
                const idType* orders
            ) const {

                // pointer used to compare against representative
                auto* currentPropagation = &propagation;

                // frequently used propagation members
                const idType& extremumIndex = currentPropagation->extremumIndex;
                auto* queue = &currentPropagation->queue;

                // add extremumIndex to queue
                queue->emplace(orders[extremumIndex],extremumIndex);
                queueMask[extremumIndex] = extremumIndex;

                // grow region until it reaches a saddle and then decide if it should continue
                while(!queue->empty()){
                    const idType v = std::get<1>(queue->top());
                    queue->pop();

                    // continue if this thread has already seen this vertex
                    if(propagationMask[v]!=nullptr)
                        continue;

                    // add neighbors to queue AND check if v is a saddle
                    bool isSaddle = false;
                    idType nNeighbors = triangulation->getVertexNeighborNumber( v );

                    idType numberOfLargerNeighbors = 0;
                    idType numberOfLargerNeighborsThisThreadVisited = 0;
                    for(idType n=0; n<nNeighbors; n++){
                        idType u;
                        triangulation->getVertexNeighbor(v,n,u);

                        // if larger neighbor
                        if( orders[u]>orders[v] ){
                            numberOfLargerNeighbors++;

                            if(propagationMask[u]==nullptr || currentPropagation!=propagationMask[u]->find())
                                isSaddle = true;
                            else
                                numberOfLargerNeighborsThisThreadVisited++;

                        } else if(queueMask[u] != extremumIndex){
                            queue->emplace(orders[u],u);
                            queueMask[u] = extremumIndex;
                        }
                    }

                    // if v is a saddle we have to check if the current thread is the last visitor
                    if(isSaddle){
                        currentPropagation->lastEncounteredCriticalPoint = v;
                        currentPropagation->terminated = 1;

                        idType numberOfRegisteredLargerVertices=0;
                        #pragma omp atomic capture
                        {
                            saddleMask[v] -= numberOfLargerNeighborsThisThreadVisited;
                            numberOfRegisteredLargerVertices = saddleMask[v];
                        }

                        // if this thread did not register the last remaining larger vertices then terminate propagation
                        if(numberOfRegisteredLargerVertices != -numberOfLargerNeighbors-1)
                            return 1;

                        // Otherwise merge propagation data
                        for(idType n=0; n<nNeighbors; n++){
                            idType u;
                            triangulation->getVertexNeighbor(v,n,u);
                            if(orders[v]<orders[u] && currentPropagation!=propagationMask[u]->find()){
                                currentPropagation = Propagation<idType>::unify(
                                    currentPropagation,
                                    propagationMask[u]
                                );
                                queue = &currentPropagation->queue;
                            }
                        }
                    }

                    // mark vertex as visited and continue
                    propagationMask[v] = currentPropagation;
                    currentPropagation->regionSize++;
                }

                this->printErr("propagation reached global minimum");

                return 0;
            }

            template<typename idType>
            int computePropagationII(
                idType* saddleMask, // used here to store registered larger vertices
                idType* queueMask, // used to mark vertices that have already been added to the queue by this thread
                Propagation<idType>** propagationMask,
                Propagation<idType>& propagation,
                idType& currentStatus,

                const ttk::Triangulation* triangulation,
                const idType* orders
            ) const {

                // pointer used to compare against representative
                auto* currentPropagation = &propagation;

                // frequently used propagation members
                const idType& extremumIndex = currentPropagation->extremumIndex;
                auto* queue = &currentPropagation->queue;

                // add extremumIndex to queue
                queue->emplace(orders[extremumIndex],extremumIndex);
                queueMask[extremumIndex] = extremumIndex;

                idType interval =9999999;

                // grow region until it reaches a saddle and then decide if it should continue
                while(!queue->empty()){
                    const idType v = std::get<1>(queue->top());
                    queue->pop();

                    // continue if this thread has already seen this vertex
                    if(propagationMask[v]!=nullptr)
                        continue;

                    const idType& orderV = orders[v];

                    // add neighbors to queue AND check if v is a saddle
                    bool isSaddle = false;
                    idType nNeighbors = triangulation->getVertexNeighborNumber( v );

                    idType numberOfLargerNeighbors = 0;
                    idType numberOfLargerNeighborsThisThreadVisited = 0;
                    for(idType n=0; n<nNeighbors; n++){
                        idType u;
                        triangulation->getVertexNeighbor(v,n,u);

                        const idType& orderU = orders[u];

                        // if larger neighbor
                        if( orderU>orderV ){
                            numberOfLargerNeighbors++;

                            if(propagationMask[u]==nullptr || currentPropagation!=propagationMask[u]->find())
                                isSaddle = true;
                            else
                                numberOfLargerNeighborsThisThreadVisited++;

                        } else if(queueMask[u] != extremumIndex){
                            queue->emplace(orderU,u);
                            queueMask[u] = extremumIndex;
                        }
                    }

                    // if v is a saddle we have to check if the current thread is the last visitor
                    if(isSaddle){
                        currentPropagation->lastEncounteredCriticalPoint = v;

                        idType numberOfRegisteredLargerVertices=0;
                        #pragma omp atomic capture
                        {
                            saddleMask[v] -= numberOfLargerNeighborsThisThreadVisited;
                            numberOfRegisteredLargerVertices = saddleMask[v];
                        }

                        // if this thread did not register the last remaining larger vertices then terminate propagation
                        if(numberOfRegisteredLargerVertices != -numberOfLargerNeighbors-1){
                            currentPropagation->terminated = 1;
                            return 1;
                        }

                        // Otherwise merge propagation data
                        for(idType n=0; n<nNeighbors; n++){
                            idType u;
                            triangulation->getVertexNeighbor(v,n,u);

                            if(orderV<orders[u] && currentPropagation!=propagationMask[u]->find()){
                                currentPropagation = Propagation<idType>::unify(
                                    currentPropagation,
                                    propagationMask[u]
                                );
                                queue = &currentPropagation->queue;
                            }
                        }
                    }

                    // mark vertex as visited and continue
                    propagationMask[v] = currentPropagation;
                    currentPropagation->regionSize++;

                    if(interval++>1000){
                        #pragma omp atomic write
                        currentStatus = orderV;
                        interval = 0;
                    }
                }

                this->printErr("propagation reached global minimum");

                return 0;
            }

            template<typename idType>
            int computePropagationIII(
                idType* saddleMask, // used here to store registered larger vertices
                idType* queueMask, // used to mark vertices that have already been added to the queue by this thread
                Propagation<idType>** propagationMask,
                Propagation<idType>& propagation,
                idType& nActivePropagations,

                const ttk::Triangulation* triangulation,
                const idType* orders,
                const idType& escapeInterval
            ) const {

                // pointer used to compare against representative
                auto* currentPropagation = &propagation;

                // frequently used propagation members
                const idType& extremumIndex = currentPropagation->extremumIndex;
                currentPropagation->lastEncounteredCriticalPoint = extremumIndex;
                auto* queue = &currentPropagation->queue;

                // a vector that will hold saddle propagations
                std::vector<Propagation<idType>*> saddlePropagations(32,nullptr);

                // add extremumIndex to queue
                queue->emplace(orders[extremumIndex],extremumIndex);
                queueMask[extremumIndex] = extremumIndex;

                idType counter = 0;

                // grow region until it reaches a saddle and then decide if it should continue
                idType v = -1;
                while(!queue->empty()){
                    v = std::get<1>(queue->top());
                    queue->pop();

                    // continue if this thread has already seen this vertex
                    if(propagationMask[v]!=nullptr)
                        continue;

                    const idType& orderV = orders[v];

                    // add neighbors to queue AND check if v is a saddle
                    bool isSaddle = false;
                    const idType nNeighbors = triangulation->getVertexNeighborNumber( v );

                    idType numberOfLargerNeighbors = 0;
                    idType numberOfLargerNeighborsThisThreadVisited = 0;
                    for(idType n=0; n<nNeighbors; n++){
                        idType u;
                        triangulation->getVertexNeighbor(v,n,u);

                        const idType& orderU = orders[u];

                        // if larger neighbor
                        if( orderU>orderV ){
                            numberOfLargerNeighbors++;

                            if(propagationMask[u]==nullptr || currentPropagation!=propagationMask[u]->find())
                                isSaddle = true;
                            else
                                numberOfLargerNeighborsThisThreadVisited++;
                        }
                        else if(queueMask[u] != extremumIndex) {
                            queue->emplace(orderU,u);
                            queueMask[u] = extremumIndex;
                        }
                    }

                    // if v is a saddle we have to check if the current thread is the last visitor
                    if(isSaddle){
                        currentPropagation->lastEncounteredCriticalPoint = v;
                        currentPropagation->terminated = 1;

                        idType numberOfRegisteredLargerVertices=0;
                        #pragma omp atomic capture
                        {
                            saddleMask[v] -= numberOfLargerNeighborsThisThreadVisited;
                            numberOfRegisteredLargerVertices = saddleMask[v];
                        }

                        // if this thread did not register the last remaining larger vertices then terminate propagation
                        if(numberOfRegisteredLargerVertices != -numberOfLargerNeighbors-1)
                            return 1;

                        // get most persistent branch
                        this->getSaddlePropagations<idType>(
                            saddlePropagations,
                            currentPropagation,
                            propagationMask,

                            triangulation,
                            v,
                            orderV,
                            nNeighbors,
                            orders
                        );

                        // merge other branches into most persistent branch
                        this->mergeSaddlePropagations<idType>(
                            saddlePropagations,
                            currentPropagation
                        );

                        // get correct queue
                        queue = &currentPropagation->queue;
                    }

                    // mark vertex as visited and continue
                    currentPropagation->regionSize++;
                    propagationMask[v] = currentPropagation;

                    if(counter++>escapeInterval){
                        counter = 0;

                        idType nActivePropagations_;
                        #pragma omp atomic read
                        nActivePropagations_ = nActivePropagations;

                        if(nActivePropagations_==1){
                            currentPropagation->terminated = 1;
                            return 1;
                        }
                    }
                }

                // if thread reached the global minimum finish propagation
                currentPropagation->terminated = 1;
                currentPropagation->lastEncounteredCriticalPoint = v;

                return 1;
            }

            // TODO: NEEDS LOGIC UPDATE
            template<typename idType, typename dataType>
            int computePropagationIV(
                idType* saddleMask,
                idType* queueMask,
                Propagation<idType>** propagationMask,
                Propagation<idType>& propagation,
                idType& nActivePropagations,
                idType* regionMask,
                idType* localOrders,
                idType* distanceField,

                const ttk::Triangulation* triangulation,
                const idType* inputOrders,
                const dataType* scalars,
                const dataType& persistenceThreshold,
                const bool& useRegionBasedIterations
            ) const {

                return 0;

                // pointer used to compare against representative
                auto* currentPropagation = &propagation;

                // frequently used propagation members
                const idType& extremumIndex = currentPropagation->extremumIndex;
                currentPropagation->lastEncounteredCriticalPoint = extremumIndex;
                auto* queue = &currentPropagation->queue;

                // largest maximum
                dataType elderScalar = scalars[extremumIndex];

                // a vector that will hold saddle propagations
                std::vector<Propagation<idType>*> saddlePropagations(32,nullptr);

                // add extremumIndex to queue
                queue->emplace(inputOrders[extremumIndex],extremumIndex);
                queueMask[extremumIndex] = extremumIndex;

                idType counter = 0;

                // grow region until it reaches a saddle and then decide if it should continue
                idType v = -1;
                while(!queue->empty()){
                    v = std::get<1>(queue->top());
                    queue->pop();

                    // continue if this thread has already seen this vertex
                    if(propagationMask[v]!=nullptr)
                        continue;

                    const idType& orderV = inputOrders[v];

                    // add neighbors to queue AND check if v is a saddle
                    bool isSaddle = false;
                    const idType nNeighbors = triangulation->getVertexNeighborNumber( v );

                    idType numberOfLargerNeighbors = 0;
                    idType numberOfLargerNeighborsThisThreadVisited = 0;
                    for(idType n=0; n<nNeighbors; n++){
                        idType u;
                        triangulation->getVertexNeighbor(v,n,u);

                        const idType& orderU = inputOrders[u];

                        // if larger neighbor
                        if( orderU>orderV ){
                            numberOfLargerNeighbors++;

                            if(propagationMask[u]==nullptr || currentPropagation!=propagationMask[u]->find())
                                isSaddle = true;
                            else
                                numberOfLargerNeighborsThisThreadVisited++;
                        }
                        else if(queueMask[u] != extremumIndex) {
                            queue->emplace(orderU,u);
                            queueMask[u] = extremumIndex;
                        }
                    }

                    // if v is a saddle we have to check if the current thread is the last visitor
                    if(isSaddle){
                        if(currentPropagation->persistent==0){
                            const dataType persistence = elderScalar>scalars[v]
                                ? elderScalar-scalars[v]
                                : scalars[v]-elderScalar;
                            if(persistence>persistenceThreshold){
                                currentPropagation->persistent = 1;

                                for(auto* c : currentPropagation->childBranches){
                                    // skip if childBranch is persistent
                                    if(c->persistent==1 || c->simplified==1)
                                        continue;

                                    // mark as simplified
                                    c->simplified = 1;

                                    // change uf for region computation
                                    c->setParentRecursive(c);

                                    // compute region
                                    this->computeRegion<idType>(
                                        regionMask,
                                        propagationMask,
                                        c,

                                        triangulation
                                    );

                                    // undo union find change
                                    c->setParentRecursive(currentPropagation);

                                    // spawn task for simplification
                                    #pragma omp task untied firstprivate(c)
                                    this->computeLocalOrderOfRegion<idType>(
                                        localOrders,
                                        distanceField,

                                        c,
                                        triangulation,
                                        regionMask,
                                        inputOrders,
                                        useRegionBasedIterations
                                    );

                                }
                            }
                        }

                        currentPropagation->lastEncounteredCriticalPoint = v;
                        currentPropagation->terminated = 1;

                        idType numberOfRegisteredLargerVertices=0;
                        #pragma omp atomic capture
                        {
                            saddleMask[v] -= numberOfLargerNeighborsThisThreadVisited;
                            numberOfRegisteredLargerVertices = saddleMask[v];
                        }

                        // if this thread did not register the last remaining larger vertices then terminate propagation
                        if(numberOfRegisteredLargerVertices != -numberOfLargerNeighbors-1){
                            #pragma omp atomic update
                            nActivePropagations--;

                            return 1;
                        }

                        // get most persistent branch
                        this->getSaddlePropagations<idType>(
                            saddlePropagations,
                            currentPropagation,
                            propagationMask,

                            triangulation,
                            v,
                            orderV,
                            nNeighbors,
                            inputOrders
                        );

                        // if the current propagation is persistent then simplify other branches
                        if(currentPropagation->persistent==1){
                            for(auto* saddlePropagation : saddlePropagations){
                                if(saddlePropagation && saddlePropagation->persistent==0 && saddlePropagation->simplified==0){
                                    saddlePropagation->simplified = 1;

                                    this->computeRegion<idType>(
                                        regionMask,
                                        propagationMask,
                                        saddlePropagation,

                                        triangulation
                                    );

                                    #pragma omp task untied firstprivate(saddlePropagation)
                                    this->computeLocalOrderOfRegion<idType>(
                                        localOrders,
                                        distanceField,

                                        saddlePropagation,
                                        triangulation,
                                        regionMask,
                                        inputOrders,
                                        useRegionBasedIterations
                                    );
                                }
                            }
                        }


                        // merge other branches into most persistent branch
                        this->mergeSaddlePropagations<idType>(
                            saddlePropagations,
                            currentPropagation
                        );

                        queue = &currentPropagation->queue;
                        elderScalar = scalars[currentPropagation->extremumIndex];
                    }

                    // mark vertex as visited and continue
                    currentPropagation->regionSize++;
                    propagationMask[v] = currentPropagation;

                    if(counter++>1000){
                        counter = 0;

                        idType nActivePropagations_;
                        #pragma omp atomic read
                        nActivePropagations_ = nActivePropagations;

                        if(nActivePropagations_==1){
                            currentPropagation->terminated = 1;
                            return 1;
                        }
                    }
                }

                // if thread reached the global minimum finish propagation
                currentPropagation->terminated = 1;
                currentPropagation->lastEncounteredCriticalPoint = v;

                return 1;
            }

            template<typename idType, typename dataType>
            int computePropagationV(
                idType* saddleMask,
                idType* queueMask,
                Propagation<idType>** propagationMask,
                Propagation<idType>& propagation,

                const ttk::Triangulation* triangulation,
                const idType* inputOrders,
                const dataType* scalars,
                const dataType& persistenceThreshold,
                const idType& escapeInterval
            ) const {

                // pointer used to compare against representative
                auto* currentPropagation = &propagation;

                // frequently used propagation members
                const idType& extremumIndex = currentPropagation->extremumIndex;
                currentPropagation->lastEncounteredCriticalPoint = extremumIndex;
                auto* queue = &currentPropagation->queue;

                // largest maximum
                dataType elderScalar = scalars[extremumIndex];

                // a vector that will hold saddle propagations
                std::vector<Propagation<idType>*> saddlePropagations(32,nullptr);

                // add extremumIndex to queue
                queue->emplace(inputOrders[extremumIndex],extremumIndex);
                queueMask[extremumIndex] = extremumIndex;

                idType counter = 0;

                // grow region until it reaches a saddle and then decide if it should continue
                idType v = -1;
                while(!queue->empty()){
                    v = std::get<1>(queue->top());
                    queue->pop();

                    // continue if this thread has already seen this vertex
                    if(propagationMask[v]!=nullptr)
                        continue;

                    const idType& orderV = inputOrders[v];

                    // add neighbors to queue AND check if v is a saddle
                    bool isSaddle = false;
                    const idType nNeighbors = triangulation->getVertexNeighborNumber( v );

                    idType numberOfLargerNeighbors = 0;
                    idType numberOfLargerNeighborsThisThreadVisited = 0;
                    for(idType n=0; n<nNeighbors; n++){
                        idType u;
                        triangulation->getVertexNeighbor(v,n,u);

                        const idType& orderU = inputOrders[u];

                        // if larger neighbor
                        if( orderU>orderV ){
                            numberOfLargerNeighbors++;

                            if(propagationMask[u]==nullptr || currentPropagation!=propagationMask[u]->find())
                                isSaddle = true;
                            else
                                numberOfLargerNeighborsThisThreadVisited++;
                        }
                        else if(queueMask[u] != extremumIndex) {
                            queue->emplace(orderU,u);
                            queueMask[u] = extremumIndex;
                        }
                    }

                    // if v is a saddle we have to check if the current thread is the last visitor
                    if(isSaddle){
                        currentPropagation->lastEncounteredCriticalPoint = v;
                        currentPropagation->terminated = 1;

                        const dataType persistence = elderScalar>scalars[v]
                            ? elderScalar-scalars[v]
                            : scalars[v]-elderScalar;
                        currentPropagation->persistent = persistence>persistenceThreshold ? 1 : 0;

                        idType numberOfRegisteredLargerVertices=0;
                        #pragma omp atomic capture
                        {
                            saddleMask[v] -= numberOfLargerNeighborsThisThreadVisited;
                            numberOfRegisteredLargerVertices = saddleMask[v];
                        }

                        // if this thread did not register the last remaining larger vertices then terminate propagation
                        if( numberOfRegisteredLargerVertices != -numberOfLargerNeighbors-1 )
                            return 1;

                        // get most persistent branch
                        this->getSaddlePropagations<idType>(
                            saddlePropagations,
                            currentPropagation,
                            propagationMask,

                            triangulation,
                            v,
                            orderV,
                            nNeighbors,
                            inputOrders
                        );

                        // merge other branches into most persistent branch
                        this->mergeSaddlePropagations<idType>(
                            saddlePropagations,
                            currentPropagation
                        );

                        queue = &currentPropagation->queue;
                        elderScalar = scalars[currentPropagation->extremumIndex];
                    }

                    // mark vertex as visited and continue
                    currentPropagation->regionSize++;
                    propagationMask[v] = currentPropagation;

                    if(counter++>escapeInterval){
                        counter = 0;
                        const dataType persistence = elderScalar>scalars[v]
                            ? elderScalar-scalars[v]
                            : scalars[v]-elderScalar;
                        currentPropagation->persistent = persistence>persistenceThreshold ? 1 : 0;
                    }

                    if(currentPropagation->persistent==1)
                        return 1;
                }

                // if thread reached the global minimum finish propagation
                currentPropagation->terminated = 1;
                currentPropagation->lastEncounteredCriticalPoint = v;

                return 1;
            }

            // propagate until trunk mode, or when getting persistent spawn new task
            template<typename idType, typename dataType>
            int computePropagationVI(
                idType* saddleMask,
                idType* queueMask,
                Propagation<idType>** propagationMask,
                Propagation<idType>& propagation,
                Propagation<idType>*& earlyEscapedPropagation,

                const ttk::Triangulation* triangulation,
                const idType* orders,
                const dataType* scalars,
                const dataType& persistenceThreshold,
                const idType& escapeInterval,
                const idType& nActivePropagations
            ) const {

                // pointer used to compare against representative
                auto* currentPropagation = &propagation;

                // frequently used propagation members
                const idType& extremumIndex = currentPropagation->extremumIndex;
                currentPropagation->lastEncounteredCriticalPoint = extremumIndex;
                auto* queue = &currentPropagation->queue;

                // largest maximum
                dataType elderScalar = scalars[extremumIndex];

                // a vector that will hold saddle propagations
                std::vector<Propagation<idType>*> saddlePropagations(32,nullptr);

                // add extremumIndex to queue
                queue->emplace(orders[extremumIndex],extremumIndex);
                queueMask[extremumIndex] = extremumIndex;

                idType counter = 0;

                // grow region until it reaches a saddle and then decide if it should continue
                idType v = -1;
                while(!queue->empty()){
                    v = std::get<1>(queue->top());
                    queue->pop();

                    // continue if this thread has already seen this vertex
                    if(propagationMask[v]!=nullptr)
                        continue;

                    const idType& orderV = orders[v];

                    // add neighbors to queue AND check if v is a saddle
                    bool isSaddle = false;
                    const idType nNeighbors = triangulation->getVertexNeighborNumber( v );

                    idType numberOfLargerNeighbors = 0;
                    idType numberOfLargerNeighborsThisThreadVisited = 0;
                    for(idType n=0; n<nNeighbors; n++){
                        idType u;
                        triangulation->getVertexNeighbor(v,n,u);

                        const idType& orderU = orders[u];

                        // if larger neighbor
                        if( orderU>orderV ){
                            numberOfLargerNeighbors++;

                            if(propagationMask[u]==nullptr || currentPropagation!=propagationMask[u]->find())
                                isSaddle = true;
                            else
                                numberOfLargerNeighborsThisThreadVisited++;
                        }
                        else if(queueMask[u] != extremumIndex) {
                            queue->emplace(orderU,u);
                            queueMask[u] = extremumIndex;
                        }
                    }

                    // if v is a saddle check if the current thread is the last visitor
                    if(isSaddle){
                        currentPropagation->lastEncounteredCriticalPoint = v;
                        currentPropagation->terminated = 1;

                        const dataType persistence = elderScalar>scalars[v]
                            ? elderScalar-scalars[v]
                            : scalars[v]-elderScalar;
                        currentPropagation->persistent = persistence>persistenceThreshold ? 1 : 0;

                        idType numberOfRegisteredLargerVertices=0;
                        #pragma omp atomic capture
                        {
                            saddleMask[v] -= numberOfLargerNeighborsThisThreadVisited;
                            numberOfRegisteredLargerVertices = saddleMask[v];
                        }

                        // if this thread did not register the last remaining larger vertices then terminate propagation
                        if( numberOfRegisteredLargerVertices != -numberOfLargerNeighbors-1 )
                            return 1;

                        // get most persistent branch
                        this->getSaddlePropagations<idType>(
                            saddlePropagations,
                            currentPropagation,
                            propagationMask,

                            triangulation,
                            v,
                            orderV,
                            nNeighbors,
                            orders
                        );

                        // merge other branches into most persistent branch
                        this->mergeSaddlePropagations<idType>(
                            saddlePropagations,
                            currentPropagation
                        );

                        queue = &currentPropagation->queue;
                        elderScalar = scalars[currentPropagation->extremumIndex];
                    }

                    // mark vertex as visited and continue
                    currentPropagation->regionSize++;
                    propagationMask[v] = currentPropagation;

                    if(counter++>escapeInterval){
                        counter = 0;

                        // check if persistence threshold is reached
                        const dataType persistence = elderScalar>scalars[v]
                            ? elderScalar-scalars[v]
                            : scalars[v]-elderScalar;
                        currentPropagation->persistent = persistence>persistenceThreshold ? 1 : 0;

                        // check if trunk mode can be started
                        idType nActivePropagations_ = 0;
                        #pragma omp atomic read
                        nActivePropagations_ = nActivePropagations;

                        if(nActivePropagations_==1)
                            return 1;
                    }

                    if(currentPropagation->persistent==1){
                        earlyEscapedPropagation = currentPropagation;
                        return 1;
                    }
                }

                // force that a propagation that reaches the global minimum is always persistent
                currentPropagation->persistent = 1;

                // if thread reached the global minimum finish propagation
                currentPropagation->terminated = 1;
                currentPropagation->lastEncounteredCriticalPoint = v;

                return 1;
            }

            // continue propagations until all other propagations are persistent or we hit trunk mode
            template<typename idType>
            int computePropagationVII(
                idType* saddleMask,
                idType* queueMask,
                Propagation<idType>** propagationMask,
                Propagation<idType>& propagation,

                const ttk::Triangulation* triangulation,
                const idType* orders,
                const idType& escapeInterval,
                const idType& nFirstPhasePropagations
            ) const {

                // pointer used to compare against representative
                auto* currentPropagation = &propagation;

                // frequently used propagation members
                const idType& extremumIndex = currentPropagation->extremumIndex;
                auto* queue = &currentPropagation->queue;

                // a vector that will hold saddle propagations
                std::vector<Propagation<idType>*> saddlePropagations(32,nullptr);

                // if this is a new Propagation
                if(currentPropagation->extremumIndex<0){
                    currentPropagation->lastEncounteredCriticalPoint = extremumIndex;

                    // add extremumIndex to queue
                    queue->emplace(orders[extremumIndex],extremumIndex);
                    queueMask[extremumIndex] = extremumIndex;
                }

                idType counter = 0;

                // grow region until it reaches a saddle and then decide if it should continue
                idType v = -1;
                while(!queue->empty()){
                    v = std::get<1>(queue->top());
                    queue->pop();

                    // continue if this thread has already seen this vertex
                    if(propagationMask[v]!=nullptr)
                        continue;

                    const idType& orderV = orders[v];

                    // add neighbors to queue AND check if v is a saddle
                    bool isSaddle = false;
                    const idType nNeighbors = triangulation->getVertexNeighborNumber( v );

                    idType numberOfLargerNeighbors = 0;
                    idType numberOfLargerNeighborsThisThreadVisited = 0;
                    for(idType n=0; n<nNeighbors; n++){
                        idType u;
                        triangulation->getVertexNeighbor(v,n,u);

                        const idType& orderU = orders[u];

                        // if larger neighbor
                        if( orderU>orderV ){
                            numberOfLargerNeighbors++;

                            if(propagationMask[u]==nullptr || currentPropagation!=propagationMask[u]->find())
                                isSaddle = true;
                            else
                                numberOfLargerNeighborsThisThreadVisited++;
                        }
                        else if(queueMask[u] != extremumIndex) {
                            queue->emplace(orderU,u);
                            queueMask[u] = extremumIndex;
                        }
                    }

                    // if v is a saddle we have to check if the current thread is the last visitor
                    if(isSaddle){
                        currentPropagation->lastEncounteredCriticalPoint = v;
                        currentPropagation->terminated = 1;

                        idType numberOfRegisteredLargerVertices=0;
                        #pragma omp atomic capture
                        {
                            saddleMask[v] -= numberOfLargerNeighborsThisThreadVisited;
                            numberOfRegisteredLargerVertices = saddleMask[v];
                        }

                        // if this thread did not register the last remaining larger vertices then terminate propagation
                        if( numberOfRegisteredLargerVertices != -numberOfLargerNeighbors-1 )
                            return 1;

                        // get most persistent branch
                        this->getSaddlePropagations<idType>(
                            saddlePropagations,
                            currentPropagation,
                            propagationMask,

                            triangulation,
                            v,
                            orderV,
                            nNeighbors,
                            orders
                        );

                        // merge other branches into most persistent branch
                        this->mergeSaddlePropagations<idType>(
                            saddlePropagations,
                            currentPropagation
                        );

                        queue = &currentPropagation->queue;
                    }

                    // mark vertex as visited and continue
                    currentPropagation->regionSize++;
                    propagationMask[v] = currentPropagation;

                    // check if all other propagations reached a persistent branch
                    if(counter++>escapeInterval){
                        counter = 0;

                        idType nFirstPhasePropagations_ = 0;
                        #pragma omp atomic read
                        nFirstPhasePropagations_ = nFirstPhasePropagations;

                        if(nFirstPhasePropagations_<1)
                            return 1;
                    }
                }

                // if thread reached the global minimum finish propagation
                currentPropagation->terminated = 1;
                currentPropagation->lastEncounteredCriticalPoint = v;

                return 1;
            }

            template<typename idType>
            int computeTrunk(
                Propagation<idType>** propagationMask,
                std::vector<Propagation<idType>>& propagations,

                const ttk::Triangulation* triangulation,
                const idType* saddleMask,
                const idType* orders,
                const idType* sortedIndices
            ) const {
                ttk::Timer timer;
                this->printMsg(
                    "Computing trunk",
                    0, 0, this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                const idType nVertices = triangulation->getNumberOfVertices();
                const idType nPropagations = propagations.size();

                Propagation<idType>* currentPropagation = nullptr;
                idType trunkIndex;

                // get largest unfinished propagation
                for(trunkIndex=nVertices-1; trunkIndex>=0; trunkIndex--){
                    const idType& v = sortedIndices[trunkIndex];

                    if(propagationMask[v]!=nullptr)
                        continue;

                    idType nNeighbors = triangulation->getVertexNeighborNumber(v);
                    for(idType n=0; n<nNeighbors; n++){
                        idType u;
                        triangulation->getVertexNeighbor(v,n,u);

                        if(propagationMask[u]!=nullptr && (currentPropagation==nullptr || orders[currentPropagation->extremumIndex]<orders[propagationMask[u]->find()->extremumIndex] )){
                            currentPropagation = propagationMask[u]->find();
                        }
                    }

                    if(currentPropagation==nullptr){
                        this->printErr("WHAT");
                        return 0;
                    }
                    break;
                }

                idType nSaddles = 0;
                idType nTrunkVertices = 0;
                if(currentPropagation!=nullptr && currentPropagation->persistent==0){

                    std::vector<Propagation<idType>*> saddlePropagations(32,nullptr);

                    // continue propagation in sorted order
                    for(; trunkIndex>=0; trunkIndex--){
                        const idType& v = sortedIndices[trunkIndex];

                        if(propagationMask[v]!=nullptr)
                            continue;

                        nTrunkVertices++;

                        // if v is a saddle
                        if(saddleMask[v]<-1){
                            nSaddles++;

                            currentPropagation->lastEncounteredCriticalPoint = v;
                            currentPropagation->terminated = 1;

                            const idType nNeighbors = triangulation->getVertexNeighborNumber(v);

                            // get most persistent branch
                            this->getSaddlePropagations<idType>(
                                saddlePropagations,
                                currentPropagation,
                                propagationMask,

                                triangulation,
                                v,
                                orders[v],
                                nNeighbors,
                                orders
                            );

                            // merge other branches into most persistent branch
                            this->mergeSaddlePropagations<idType>(
                                saddlePropagations,
                                currentPropagation
                            );
                        }

                        propagationMask[v] = currentPropagation;
                        currentPropagation->regionSize++;
                    }

                    // mark main branch as terminated
                    currentPropagation->terminated = 1;
                }

                std::stringstream vFraction;
                vFraction << std::fixed << std::setprecision(2) << ((float)nTrunkVertices/(float)nVertices);

                this->printMsg(
                    "Computing trunk ("+std::to_string(nSaddles)+"|"+std::to_string(nTrunkVertices)+"|"+vFraction.str()+")",
                    1, timer.getElapsedTime(), this->threadNumber_
                );

                return 1;
            }

            template<typename idType, typename dataType>
            int computeTrunkII(
                Propagation<idType>** propagationMask,
                std::vector<Propagation<idType>>& propagations,

                const ttk::Triangulation* triangulation,
                const idType* saddleMask,
                const idType* orders,
                const std::vector<std::tuple<idType,idType,idType>>& sortedIndices,
                const dataType* scalars,
                const dataType& persistenceThreshold
            ) const {
                ttk::Timer timer;
                this->printMsg(
                    "Computing trunk",
                    0, 0, this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                const idType nVertices = triangulation->getNumberOfVertices();

                Propagation<idType>* currentPropagation = nullptr;
                idType trunkIndex;

                // get largest unfinished propagation
                for(trunkIndex=0; trunkIndex<nVertices; trunkIndex++){
                    const idType& v = std::get<2>(sortedIndices[trunkIndex]);

                    if(propagationMask[v]!=nullptr)
                        continue;

                    idType nNeighbors = triangulation->getVertexNeighborNumber(v);
                    for(idType n=0; n<nNeighbors; n++){
                        idType u;
                        triangulation->getVertexNeighbor(v,n,u);

                        if(propagationMask[u]!=nullptr && (currentPropagation==nullptr || orders[currentPropagation->extremumIndex]<orders[propagationMask[u]->find()->extremumIndex] )){
                            currentPropagation = propagationMask[u]->find();
                        }
                    }

                    break;
                }

                idType nSaddles = 0;
                idType nTrunkVertices = 0;
                if(currentPropagation!=nullptr){

                    dataType elderScalar = scalars[currentPropagation->extremumIndex];
                    std::vector<Propagation<idType>*> saddlePropagations(32,nullptr);

                    // continue propagation in sorted order
                    for(; trunkIndex<nVertices; trunkIndex++){
                        const idType& v = std::get<2>(sortedIndices[trunkIndex]);

                        if(propagationMask[v]!=nullptr)
                            continue;

                        nTrunkVertices++;

                        // if v is a saddle
                        if(saddleMask[v]<-1){
                            nSaddles++;

                            currentPropagation->lastEncounteredCriticalPoint = v;
                            currentPropagation->terminated = 1;

                            const idType nNeighbors = triangulation->getVertexNeighborNumber(v);

                            // get most persistent branch
                            this->getSaddlePropagations<idType>(
                                saddlePropagations,
                                currentPropagation,
                                propagationMask,

                                triangulation,
                                v,
                                orders[v],
                                nNeighbors,
                                orders
                            );

                            // merge other branches into most persistent branch
                            this->mergeSaddlePropagations<idType>(
                                saddlePropagations,
                                currentPropagation
                            );

                            elderScalar = scalars[currentPropagation->extremumIndex];
                        }

                        propagationMask[v] = currentPropagation;
                        currentPropagation->regionSize++;

                        const dataType persistence = elderScalar>scalars[v]
                            ? elderScalar-scalars[v]
                            : scalars[v]-elderScalar;
                        if(persistence>persistenceThreshold){
                            currentPropagation->persistent = 1;
                            break;
                        }
                    }

                    // mark main branch as terminated
                    currentPropagation->terminated = 1;
                }

                std::stringstream vFraction;
                vFraction << std::fixed << std::setprecision(2) << ((float)nTrunkVertices/(float)nVertices);

                this->printMsg(
                    "Computing trunk ("+std::to_string(nSaddles)+"|"+std::to_string(nTrunkVertices)+"|"+vFraction.str()+")",
                    1, timer.getElapsedTime(), this->threadNumber_
                );

                return 1;
            }

            template<typename idType>
            int spawnInterleavedTasks(
                idType* localOrders,
                idType* regionMask,
                Propagation<idType>** propagationMask,
                std::vector<Propagation<idType>>& propagations,
                idType* distanceField,

                const ttk::Triangulation* triangulation,
                const idType* inputOrders,
                const bool& useRegionBasedIterations,
                const std::vector<idType>& progress
            )const {
                // const idType nPropagations = propagations.size();
                // const idType nVertices = triangulation->getNumberOfVertices();

                // idType pivot = 0;
                // // idType lastActivePropagation = 0;
                // // idType sanity = 0;
                // // for(idType i=0; i<this->threadNumber_; i++){
                // //     #pragma omp atomic read
                // //     pivot=progress[i];

                // //     if(pivot!=-nVertices-1){
                // //         lastActivePropagation = i;
                // //         sanity++;
                // //     }
                // // }

                // do {
                //     #pragma omp atomic read
                //     pivot = progress[0];

                //     // get progress of slowest propagation
                //     for(idType t=1; t<this->threadNumber_; t++){
                //         idType progress_ = 0;
                //         #pragma omp atomic read
                //         progress_ = progress[t];
                //         if(pivot<progress_)
                //             pivot=progress_;
                //     }

                //     // #pragma omp atomic read
                //     // pivot = progress[lastActivePropagation];

                //     for(idType p=0; p<nPropagations; p++){
                //         auto propagation = &propagations[p];

                //         signed char terminated;
                //         // #pragma omp atomic read
                //         terminated = propagation->terminated;

                //         idType saddleIndex;
                //         // #pragma omp atomic read
                //         saddleIndex = propagation->lastEncounteredCriticalPoint;

                //         // Propagation<idType>* parent;
                //         // // #pragma omp atomic read
                //         // parent = propagation->parent;

                //         // if the propagation is complete
                //         if(terminated==0 || propagation->temp==1)
                //             continue;

                //         if(inputOrders[saddleIndex]<=pivot)
                //             continue;

                //         propagation->temp = 1;

                //         #pragma omp task firstprivate(propagation) priority(1)
                //         {
                //             this->computeRegion<idType>(
                //                 regionMask,
                //                 propagationMask,
                //                 propagation,

                //                 triangulation
                //             );

                //             this->computeLocalOrderOfRegion<idType>(
                //                 localOrders,
                //                 distanceField,

                //                 propagation,
                //                 triangulation,
                //                 regionMask,
                //                 inputOrders,
                //                 useRegionBasedIterations
                //             );
                //         }
                //     }
                // } while(pivot!=-nVertices-1);

                return 0;
            }

            template<typename idType>
            int computeInterleaving(
                idType* localOrders,
                idType* regionMask,
                idType* queueMask,
                Propagation<idType>** propagationMask,
                std::vector<Propagation<idType>>& propagations,
                idType* distanceField,

                const std::vector<idType>& unauthorizedMaxima,
                const ttk::Triangulation* triangulation,
                const idType* inputOrders,
                const bool& useRegionBasedIterations
            ) const {
                const idType nVertices = triangulation->getNumberOfVertices();
                const idType nPropagations = unauthorizedMaxima.size();

                ttk::Timer timer;
                this->printMsg(
                    "Interleaved Computation ("+std::to_string(nPropagations)+")",
                    0, 0, this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                std::vector<idType> progress(this->threadNumber_,nVertices+1);
                idType nActivePropagations = nPropagations;

                int status = 1;
                #pragma omp parallel num_threads(this->threadNumber_)
                #pragma omp single
                for(idType t=0; t<this->threadNumber_; t++){

                    #pragma omp task firstprivate(t) priority(2)
                    {
                        while(true) {
                            idType nActivePropagations_ = 0;
                            #pragma omp atomic capture
                            {
                                nActivePropagations--;
                                nActivePropagations_ = nActivePropagations;
                            }

                            if(nActivePropagations_<0){
                                progress[t] = -nVertices-1;

                                if(
                                    // (this->threadNumber_>1 && nActivePropagations_==-this->threadNumber_+1)
                                    (this->threadNumber_>1 && nActivePropagations_==-1)
                                    ||
                                    (this->threadNumber_<2)
                                ){
                                    this->spawnInterleavedTasks<idType>(
                                        localOrders,
                                        regionMask,
                                        propagationMask,
                                        propagations,
                                        distanceField,

                                        triangulation,
                                        inputOrders,
                                        useRegionBasedIterations,
                                        progress
                                    );
                                }

                                break;
                            }

                            this->computePropagationII<idType>(
                                regionMask,
                                queueMask,
                                propagationMask,
                                propagations[nActivePropagations_],
                                progress[t],

                                triangulation,
                                inputOrders
                            );
                        }
                    }
                }
                if(!status)
                    return 0;

                this->printMsg(
                    "Interleaved Computation ("+std::to_string(nPropagations)+")",
                    1, timer.getElapsedTime(), this->threadNumber_
                );

                return 1;
            }

            template<typename idType>
            int sortPropagations(
                std::vector<Propagation<idType>>& propagations,

                const idType* orders
            ) const {
                ttk::Timer timer;

                const idType nPropagations = propagations.size();

                this->printMsg(
                    "Sort propagations ("+std::to_string(nPropagations)+")",
                    0, timer.getElapsedTime(), this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                std::vector<std::tuple<idType,idType>> sortedPropagations(nPropagations);
                #pragma omp parallel for num_threads(this->threadNumber_)
                for(idType p=0; p<nPropagations; p++){
                    std::get<0>(sortedPropagations[p]) = -orders[propagations[p].extremumIndex];
                    std::get<1>(sortedPropagations[p]) = propagations[p].extremumIndex;
                }

                {
                    #ifdef TTK_ENABLE_OPENMP
                        #ifdef __clang__
                            this->printWrn("Caution, outside GCC, sequential sort");
                            std::sort(sortedPropagations.begin(), sortedPropagations.end());
                        #else
                            omp_set_num_threads(this->threadNumber_);
                            __gnu_parallel::sort(sortedPropagations.begin(), sortedPropagations.end());
                            omp_set_num_threads(1);
                        #endif
                    #else
                        this->printWrn("Caution, outside GCC, sequential sort");
                        std::sort(sortedPropagations.begin(), sortedPropagations.end());
                    #endif
                }

                #pragma omp parallel for num_threads(this->threadNumber_)
                for(idType p=0; p<nPropagations; p++){
                    propagations[p].extremumIndex = std::get<1>(sortedPropagations[p]);
                }

                this->printMsg(
                    "Sort propagations ("+std::to_string(nPropagations)+")",
                    1, timer.getElapsedTime(), this->threadNumber_
                );

                return 1;
            };

            template<typename idType, typename dataType>
            int initializeScalars(
                dataType* outputScalars,

                const dataType* inputScalars,
                const idType& nVertices
            ) const {
                ttk::Timer timer;

                this->printMsg(
                    "Initialize scalars ("+std::to_string(nVertices)+")",
                    0, 0, this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                // init region/queue/propagation mask
                #pragma omp parallel for num_threads(this->threadNumber_)
                for(idType v=0; v<nVertices; v++){
                    outputScalars[v] = inputScalars[v];
                }

                this->printMsg(
                    "Initialize scalars ("+std::to_string(nVertices)+")",
                    1, timer.getElapsedTime(), this->threadNumber_
                );

                return 1;
            };

            template<typename idType>
            int initializePropagations(
                std::vector<Propagation<idType>>& propagations,
                idType* saddleOrRegionMask,
                idType* queueMask,
                Propagation<idType>** propagationMask,

                const std::vector<idType>& unauthorizedExtrema,
                const idType& nVertices,

                idType* localOrders = nullptr
            ) const {
                ttk::Timer timer;

                const idType nPropagations = unauthorizedExtrema.size();

                this->printMsg(
                    "Initialize propagations ("+std::to_string(nPropagations)+")",
                    0, timer.getElapsedTime(), this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                // init region/queue/propagation mask
                if(localOrders){
                    #pragma omp parallel for num_threads(this->threadNumber_)
                    for(idType i=0; i<nVertices; i++){
                        saddleOrRegionMask[i] = -1;
                        localOrders[i] = 1;
                        queueMask[i] = -1;
                        propagationMask[i] = nullptr;
                    }
                } else {
                    #pragma omp parallel for num_threads(this->threadNumber_)
                    for(idType i=0; i<nVertices; i++){
                        saddleOrRegionMask[i] = -1;
                        queueMask[i] = -1;
                        propagationMask[i] = nullptr;
                    }
                }

                propagations.clear();
                propagations.resize(nPropagations);
                #pragma omp parallel for num_threads(this->threadNumber_)
                for(idType i=0; i<nPropagations; i++)
                    propagations[i].extremumIndex = unauthorizedExtrema[i];

                this->printMsg(
                    "Initialize propagations ("+std::to_string(nPropagations)+")",
                    1, timer.getElapsedTime(), this->threadNumber_
                );

                return 1;
            };

            template<typename idType>
            int finalizePropagations(
                std::vector<Propagation<idType>*>& masterPropagations,
                std::vector<Propagation<idType>>& propagations,

                const idType nVertices
            ) const {
                ttk::Timer timer;

                const idType nPropagations = propagations.size();

                this->printMsg(
                    "Finalizing propagations ("+std::to_string(nPropagations)+")",
                    0, timer.getElapsedTime(), this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                idType nRegionVertices = 0;
                idType nMasterPropagations=0;
                masterPropagations.clear();
                masterPropagations.resize(nPropagations);
                for(idType p=0; p<nPropagations; p++){
                    auto* propagation = &propagations[p];
                    if(propagation->parent == propagation && propagation->terminated==1){
                        nRegionVertices = nRegionVertices + propagation->regionSize;
                        masterPropagations[nMasterPropagations++] = propagation;
                    }
                }
                masterPropagations.resize(nMasterPropagations);

                std::stringstream pFraction, vFraction;
                pFraction << std::fixed << std::setprecision(2) << ((float)nMasterPropagations/(float)nPropagations);
                vFraction << std::fixed << std::setprecision(2) << ((float)nRegionVertices/(float)nVertices);

                this->printMsg(
                    "Finalizing propagations ("+std::to_string(nMasterPropagations)+"|"+pFraction.str()+"|"+vFraction.str()+")",
                    1, timer.getElapsedTime(), this->threadNumber_
                );

                return 1;
            };

            template<typename idType, typename dataType>
            int finalizePropagationsByPersistence(
                std::vector<Propagation<idType>*>& masterPropagations,
                std::vector<Propagation<idType>>& propagations,

                const dataType* scalars,
                const idType nVertices
            ) const {
                ttk::Timer timer;

                const idType nPropagations = propagations.size();

                this->printMsg(
                    "Finalizing propagations ("+std::to_string(nPropagations)+")",
                    0, timer.getElapsedTime(), this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                idType nRegionVertices = 0;
                idType nMasterPropagations=0;
                masterPropagations.clear();
                masterPropagations.resize(nPropagations);
                #pragma omp parallel for num_threads(threadNumber_) reduction(+:nRegionVertices)
                for(idType p=0; p<nPropagations; p++){
                    Propagation<idType>* propagation = &propagations[p];

                    // if the propagation is persistent or a is the child of a persistent branch skip
                    if(propagation->persistent==1 || (propagation->parentBranch && propagation->parentBranch->persistent==0))
                        continue;

                    propagation->setParentRecursive(propagation);
                    propagation->simplified = 1;

                    idType nMasterPropagations_;
                    #pragma omp atomic capture
                    nMasterPropagations_ = nMasterPropagations++;

                    masterPropagations[nMasterPropagations_] = propagation;

                    nRegionVertices += propagation->regionSize;
                }
                masterPropagations.resize(nMasterPropagations);

                std::stringstream pFraction, vFraction;
                pFraction << std::fixed << std::setprecision(2) << ((float)nMasterPropagations/(float)nPropagations);
                vFraction << std::fixed << std::setprecision(2) << ((float)nRegionVertices/(float)nVertices);

                this->printMsg(
                    "Finalizing propagations ("+std::to_string(nMasterPropagations)+"|"+pFraction.str()+"|"+vFraction.str()+")",
                    1, timer.getElapsedTime(), this->threadNumber_
                );

                return 1;
            };

            template<typename idType, typename dataType>
            int finalizePropagationsByPersistenceII(
                std::vector<Propagation<idType>*>& masterPropagations,
                std::vector<Propagation<idType>>& propagations,

                const ttk::Triangulation* triangulation,
                const idType* orders,
                const dataType* scalars,
                const dataType& persistenceThreshold,
                const idType nVertices
            ) const {
                ttk::Timer timer;

                const idType nPropagations = propagations.size();

                this->printMsg(
                    "Finalizing propagations ("+std::to_string(nPropagations)+")",
                    0, timer.getElapsedTime(), this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                idType nRegionVertices = 0;
                idType nMasterPropagations=0;
                masterPropagations.clear();
                masterPropagations.resize(nPropagations);
                #pragma omp parallel for num_threads(threadNumber_) reduction(+:nRegionVertices)
                for(idType p=0; p<nPropagations; p++){
                    Propagation<idType>* propagation = &propagations[p];

                    const dataType persistence = scalars[propagation->extremumIndex]>scalars[propagation->lastEncounteredCriticalPoint]
                        ? scalars[propagation->extremumIndex]-scalars[propagation->lastEncounteredCriticalPoint]
                        : scalars[propagation->lastEncounteredCriticalPoint]-scalars[propagation->extremumIndex];

                    propagation->persistent = persistence>persistenceThreshold ? 1 : 0;

                    if(propagation->persistent==1 || !propagation->parentBranch){
                        propagation->parent = propagation;
                    } else {
                        const dataType parentPersistence =
                            scalars[propagation->parentBranch->extremumIndex]>scalars[propagation->parentBranch->lastEncounteredCriticalPoint]
                            ? scalars[propagation->parentBranch->extremumIndex]-scalars[propagation->parentBranch->lastEncounteredCriticalPoint]
                            : scalars[propagation->parentBranch->lastEncounteredCriticalPoint]-scalars[propagation->parentBranch->extremumIndex];

                        if(parentPersistence<=persistenceThreshold)
                            continue;

                        // check if a propagation is still an extremum
                        bool hasLargerNeighbor = false;
                        bool hasSmallerNeighbor = false;

                        idType nNeighbors = triangulation->getVertexNeighborNumber( propagation->extremumIndex );
                        for(idType n=0; n<nNeighbors; n++){
                            idType u;
                            triangulation->getVertexNeighbor(propagation->extremumIndex,n,u);

                            if( orders[propagation->extremumIndex]<orders[u] )
                                hasLargerNeighbor = true;
                            else
                                hasSmallerNeighbor = true;
                        }

                        if(hasLargerNeighbor && hasSmallerNeighbor){
                            continue;
                        }

                        propagation->setParentRecursive(propagation);
                        propagation->simplified = 1;

                        idType nMasterPropagations_;
                        #pragma omp atomic capture
                        nMasterPropagations_ = nMasterPropagations++;

                        masterPropagations[nMasterPropagations_] = propagation;

                        nRegionVertices += propagation->regionSize;
                    }
                }
                masterPropagations.resize(nMasterPropagations);

                std::stringstream pFraction, vFraction;
                pFraction << std::fixed << std::setprecision(2) << ((float)nMasterPropagations/(float)nPropagations);
                vFraction << std::fixed << std::setprecision(2) << ((float)nRegionVertices/(float)nVertices);

                this->printMsg(
                    "Finalizing propagations ("+std::to_string(nMasterPropagations)+"|"+pFraction.str()+"|"+vFraction.str()+")",
                    1, timer.getElapsedTime(), this->threadNumber_
                );

                return 1;
            };

            template<typename idType>
            int computePropagations(
                idType* saddleMask,
                idType* queueMask,
                Propagation<idType>** propagationMask,
                std::vector<Propagation<idType>>& propagations,

                const ttk::Triangulation* triangulation,
                const idType* orders
            ) const {
                ttk::Timer timer;

                const idType nPropagations = propagations.size();
                this->printMsg(
                    "Computing propagations ("+std::to_string(nPropagations)+")",
                    0, 0, this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                int status = 1;
                // compute propagations
                #pragma omp parallel for schedule(dynamic,1) num_threads(this->threadNumber_)
                for(idType p=0; p<nPropagations; p++){
                    int localStatus = this->computePropagation<idType>(
                        saddleMask,
                        queueMask,
                        propagationMask,
                        propagations[p],

                        triangulation,
                        orders
                    );

                    if(!localStatus)
                        status = 0;
                }

                if(!status) return 0;

                this->printMsg(
                    "Computing propagations ("+std::to_string(nPropagations)+")",
                    1, timer.getElapsedTime(), this->threadNumber_
                );

                return 1;
            }

            template<typename idType>
            int computeDynamicPropagations(
                idType* saddleMask,
                idType* queueMask,
                Propagation<idType>** propagationMask,
                std::vector<Propagation<idType>>& propagations,

                const ttk::Triangulation* triangulation,
                const idType* orders,
                const idType& escapeInterval
            ) const {
                ttk::Timer timer;

                int status = 1;
                const idType nPropagations = propagations.size();
                idType nActivePropagations = nPropagations;

                this->printMsg(
                    "Computing dynamic propagations ("+std::to_string(nPropagations)+")",
                    0, 0, this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                // compute propagations
                #pragma omp parallel for schedule(dynamic,1) num_threads(this->threadNumber_)
                for(idType p=0; p<nPropagations; p++){
                    int localStatus = this->computePropagationIII<idType>(
                        saddleMask,
                        queueMask,
                        propagationMask,
                        propagations[p],
                        nActivePropagations,

                        triangulation,
                        orders,
                        escapeInterval
                    );
                    if(!localStatus)
                        status = 0;

                    #pragma omp atomic update
                    nActivePropagations--;
                }
                if(!status) return 0;

                this->printMsg(
                    "Computing dynamic propagations ("+std::to_string(nPropagations)+")",
                    1, timer.getElapsedTime(), this->threadNumber_
                );

                return 1;
            }

            template<typename idType, typename dataType>
            int computePersistentPropagations(
                idType* saddleMask,
                idType* queueMask,
                Propagation<idType>** propagationMask,
                std::vector<Propagation<idType>>& propagations,

                const ttk::Triangulation* triangulation,
                const idType* orders,
                const dataType* scalars,
                const dataType& persistenceThreshold,
                const idType& escapeInterval
            ) const {
                ttk::Timer timer;

                int status = 1;
                const idType nPropagations = propagations.size();

                this->printMsg(
                    "Computing persistent propagations ("+std::to_string(nPropagations)+")",
                    0, 0, this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                idType propagationIndex = 0;
                idType nActivePropagations = nPropagations;
                idType nFirstPhasePropagations = nPropagations;

                // compute propagations
                #pragma omp parallel num_threads(this->threadNumber_)
                #pragma omp single
                for(idType t=0; t<this->threadNumber_; t++){
                    #pragma omp task
                    {
                        idType propagationIndex_;
                        while(true){
                            #pragma omp atomic capture
                            propagationIndex_ = propagationIndex++;

                            if(propagationIndex_>=nPropagations)
                                break;

                            Propagation<idType>* earlyEscapedPropagation = nullptr;
                            int localStatus = this->computePropagationVI<idType,dataType>(
                                saddleMask,
                                queueMask,
                                propagationMask,
                                propagations[propagationIndex_],
                                earlyEscapedPropagation,

                                triangulation,
                                orders,
                                scalars,
                                persistenceThreshold,
                                escapeInterval,
                                nActivePropagations
                            );
                            if(!localStatus)
                                status = 0;

                            #pragma omp atomic update
                            nFirstPhasePropagations--;

                            if(earlyEscapedPropagation){
                                #pragma omp task firstprivate(earlyEscapedPropagation)
                                {
                                    this->computePropagationVII<idType>(
                                        saddleMask,
                                        queueMask,
                                        propagationMask,
                                        *earlyEscapedPropagation,

                                        triangulation,
                                        orders,
                                        escapeInterval,
                                        nFirstPhasePropagations
                                    );

                                    #pragma omp atomic update
                                    nActivePropagations--;
                                }
                            } else {
                                #pragma omp atomic update
                                nActivePropagations--;
                            }
                        }
                    }
                }
                if(!status) return 0;



                idType regionSizes = 0;
                if(this->debugLevel_>3){
                    double t = timer.getElapsedTime();
                    idType nVertices = triangulation->getNumberOfVertices();
                    #pragma omp parallel for num_threads(this->threadNumber_) reduction(+:regionSizes)
                    for(idType v=0; v<nVertices; v++)
                        if(propagationMask[v]!=nullptr) regionSizes ++;

                    std::stringstream vFraction;
                    vFraction << std::fixed << std::setprecision(2) << ((float)regionSizes/(float)nVertices);
                    this->printMsg(
                        "Computing persistent propagations ("+std::to_string(nPropagations)+"|"+vFraction.str()+")",
                        1, t, this->threadNumber_
                    );
                } else
                    this->printMsg(
                        "Computing persistent propagations ("+std::to_string(nPropagations)+")",
                        1, timer.getElapsedTime(), this->threadNumber_
                    );


                return 1;
            }

            template<typename idType>
            int computeLevelSetComponents(
                idType* localOrders,

                const ttk::Triangulation* triangulation,
                const idType* regionMask,
                const idType& regionID,
                const std::vector<idType>& region,
                const std::vector<idType>& seedVertices,
                const idType* distanceField,
                const idType& localOrderSortingDirection // controls if the first visited vertex is the last or first element of the local orders
            ) const {

                // NOTE: this function uses the localOrders during computation to add every region vertex exectly once
                for(const auto& v: region)
                    localOrders[v] = regionID;

                // init priority queue
                std::priority_queue<
                    std::pair<idType,idType>,
                    std::vector<std::pair<idType,idType>>
                > queue;

                for(const auto& v: seedVertices){
                    queue.emplace( distanceField[v], v );
                    localOrders[v] = -1;
                }

                std::vector<idType> localVertexSequence(region.size());

                idType q=0;
                while(!queue.empty()){
                    idType v = std::get<1>(queue.top());
                    queue.pop();

                    localVertexSequence[q++] = v;

                    idType nNeighbors = triangulation->getVertexNeighborNumber( v );
                    for(idType n=0; n<nNeighbors; n++){
                        idType u;
                        triangulation->getVertexNeighbor(v,n,u);

                        // if u has not been popped from the queue add to queue
                        if(localOrders[u]==regionID){
                            queue.emplace( distanceField[u], u );
                            localOrders[u] = -1;
                        }
                    }
                }

                idType localOrder = -1;
                if(localOrderSortingDirection<0){
                    for(idType i=0, j=region.size(); i<j; i++)
                        localOrders[ localVertexSequence[i] ] = localOrder--;
                } else {
                    for(idType i=region.size()-1; i>=0; i--)
                        localOrders[ localVertexSequence[i] ] = localOrder--;
                }

                return 1;
            }

            template<typename idType>
            int computeLocalOrderOfRegion(
                idType* localOrders,
                idType* distanceField,

                const Propagation<idType>* propagation,
                const ttk::Triangulation* triangulation,
                const idType* regionMask,
                const idType* inputOrders,
                const bool& useRegionBasedIterations
            ) const {
                propagation->nIterations = 1;

                if(propagation->regionSize==1){
                    localOrders[ propagation->region[0] ] = -1;
                    return 1;
                }

                const idType& extremumIndex = propagation->extremumIndex;
                const idType& saddleIndex = propagation->lastEncounteredCriticalPoint;

                // init distance field
                for(const auto& v: propagation->region)
                    distanceField[v] = inputOrders[v];

                // there is always only one authorized maximum, which is next to the saddle
                std::vector<idType> saddleNeighbors;
                {
                    // get all saddle neighbors inside region
                    idType nNeighbors = triangulation->getVertexNeighborNumber( saddleIndex );
                    for(idType n=0; n<nNeighbors; n++){
                        idType u;
                        triangulation->getVertexNeighbor(saddleIndex,n,u);
                        if(regionMask[u]==extremumIndex){
                            saddleNeighbors.push_back( u );
                            distanceField[u] = std::numeric_limits<idType>::max(); // force that there is always one maximum next to the saddle
                        }
                    }
                }

                int status = 1;

                // start first iteration from saddle neighbors
                status = this->computeLevelSetComponents<idType>(
                    localOrders,

                    triangulation,
                    regionMask,
                    saddleIndex,
                    propagation->region,
                    saddleNeighbors,
                    distanceField,
                    -1
                );
                if(!status)
                    return 0;

                if(!useRegionBasedIterations)
                    return 1;

                bool containsResidualExtrema = false;

                // get all vertices on the boundary
                std::vector<idType> regionBoundary(propagation->regionSize);
                {
                    idType boundaryWriteIndex = 0;

                    for(const auto& v: propagation->region){
                        bool isOnRegionBoundary = false;
                        bool hasSmallerNeighbor = false;

                        idType nNeighbors = triangulation->getVertexNeighborNumber( v );
                        for(idType n=0; n<nNeighbors; n++){
                            idType u;
                            triangulation->getVertexNeighbor(v,n,u);

                            // if u is not inside region -> v is on region boundary
                            if(regionMask[u]!=extremumIndex){
                                isOnRegionBoundary = true;
                            } else if (localOrders[u]<localOrders[v]){
                                hasSmallerNeighbor = true;
                            }
                        }

                        if(isOnRegionBoundary)
                            regionBoundary[boundaryWriteIndex++] = v;
                        else if(!hasSmallerNeighbor)
                            containsResidualExtrema = true;
                    }

                    regionBoundary.resize(boundaryWriteIndex);
                }

                idType localOrderSortingDirection = -1;
                while(containsResidualExtrema){
                    localOrderSortingDirection*=-1;

                    // set seed vertices and init distance field
                    std::vector<idType>* seedVertices;
                    if(localOrderSortingDirection>0){
                        seedVertices = &regionBoundary;
                        for(const auto& v: propagation->region)
                            distanceField[v] = -localOrders[v];
                    } else {
                        seedVertices = &saddleNeighbors;
                        for(const auto& v: propagation->region)
                            distanceField[v] = localOrders[v];
                    }

                    status = this->computeLevelSetComponents<idType>(
                        localOrders,

                        triangulation,
                        regionMask,
                        saddleIndex,
                        propagation->region,
                        *seedVertices,
                        distanceField,
                        localOrderSortingDirection
                    );
                    if(!status)
                        return 0;

                    // check if number of extrema correpsonds to number of authorized extrema
                    idType nUnauthorizedMinima = 0;
                    idType nUnauthorizedMaxima = 0;

                    for(const auto& v: propagation->region){
                        // check if v is on the boundary and if it has no smaller neighbors inside region
                        bool hasSmallerNeighbor = false;
                        bool hasLargerNeighbor = false;
                        bool isOnRegionBoundary = false;
                        bool isNextToSaddle = false;

                        idType nNeighbors = triangulation->getVertexNeighborNumber( v );
                        for(idType n=0; n<nNeighbors; n++){
                            idType u;
                            triangulation->getVertexNeighbor(v,n,u);

                            if(regionMask[u]!=extremumIndex){
                                isOnRegionBoundary = true;

                                if(u==saddleIndex)
                                    isNextToSaddle = true;

                                continue;
                            }

                            if(localOrders[u]<localOrders[v]) {
                                hasSmallerNeighbor = true;
                            } else {
                                hasLargerNeighbor = true;
                            }
                        }

                        if(!hasLargerNeighbor && !isNextToSaddle){
                            nUnauthorizedMaxima++;
                        } else if(!hasSmallerNeighbor && !isOnRegionBoundary){
                            nUnauthorizedMinima++;
                        }
                    }

                    propagation->nIterations++;

                    containsResidualExtrema = nUnauthorizedMaxima>0 || nUnauthorizedMinima>0;
                }

                return 1;
            }

            template<typename idType>
            int computeLocalOrderOfRegions(
                idType* localOrders,
                idType* distanceField,

                const ttk::Triangulation* triangulation,
                const idType* regionMask,
                const idType* inputOrders,
                const std::vector<Propagation<idType>*>& propagations,
                const bool& useRegionBasedIterations
            ) const {
                ttk::Timer timer;

                const idType nPropagations = propagations.size();
                this->printMsg(
                    "Computing local order of regions ("+std::to_string(nPropagations)+")",
                    0, 0, this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                int status = 1;
                #pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
                for(idType p=0; p<nPropagations; p++){
                    int localStatus = this->computeLocalOrderOfRegion<idType>(
                        localOrders,
                        distanceField,

                        propagations[p],
                        triangulation,
                        regionMask,
                        inputOrders,
                        useRegionBasedIterations
                    );
                    if(!localStatus)
                        status = 0;
                }
                if(!status)
                    return 0;

                if(this->debugLevel_<4){
                    this->printMsg( "Computing local order of regions ("+std::to_string(nPropagations)+")",
                        1, timer.getElapsedTime(), this->threadNumber_
                    );
                } else {
                    idType min = propagations[0]->nIterations;
                    idType max = min;
                    idType avg = 0;

                    for(idType p=0; p<nPropagations; p++){
                        const auto propagation = propagations[p];
                        if(min>propagation->nIterations)
                            min=propagation->nIterations;
                        if(max<propagation->nIterations)
                            max=propagation->nIterations;
                        avg += propagation->nIterations;
                    }

                    avg /= nPropagations;

                    this->printMsg(
                        "Computing local order of regions ("+std::to_string(nPropagations)
                            + "|" + std::to_string(min)
                            + "|" + std::to_string(avg)
                            + "|" + std::to_string(max)
                        +")",
                        1, timer.getElapsedTime(), this->threadNumber_
                    );
                }

                return 1;
            }

            template<typename idType>
            int detectAndRemoveUnauthorizedMaxima(
                std::vector<idType>& unauthorizedMaxima,
                idType* outputOrder,
                idType* localOrders,
                idType* regionMask,
                idType* queueMask,
                Propagation<idType>** propagationMask,
                std::vector<Propagation<idType>>& propagations,
                std::vector<Propagation<idType>*>& masterPropagations,
                std::vector<std::tuple<idType,idType,idType>>& sortedIndices,
                idType& nRemovedMaxima,

                const ttk::Triangulation* triangulation,
                const idType* authorizedExtremaIndices,
                const idType& nAuthorizedExtremaIndices,
                const idType* inputOrders,
                const bool&   useRegionBasedIterations,
                const bool&   useInterleaving
            ) const {
                const idType nVertices = triangulation->getNumberOfVertices();

                int status = 0;

                // Classify Critical Points
                status = this->detectUnauthorizedMaxima<idType>(
                    unauthorizedMaxima,
                    localOrders, // used here to temporarily store preservation mask

                    triangulation,
                    inputOrders,
                    authorizedExtremaIndices,
                    nAuthorizedExtremaIndices
                );
                if(!status) return 0;

                // if nothing to remove return
                if(unauthorizedMaxima.size()<1)
                    return 1;

                // init propagations
                status = this->initializePropagations<idType>(
                    propagations,
                    regionMask,
                    queueMask,
                    propagationMask,

                    unauthorizedMaxima,
                    nVertices,

                    localOrders
                );
                if(!status) return 0;

                if(!useInterleaving){
                    // compute propagations
                    status = this->computePropagations<idType>(
                        regionMask,
                        queueMask,
                        propagationMask,
                        propagations,

                        triangulation,
                        inputOrders
                    );
                    if(!status) return 0;

                    // finalize master propagations
                    status = this->finalizePropagations<idType>(
                        masterPropagations,
                        propagations,

                        nVertices
                    );
                    if(!status) return 0;
                    nRemovedMaxima = masterPropagations.size();

                    // compute regions
                    status = this->computeRegions<idType>(
                        regionMask,
                        propagationMask,
                        masterPropagations,

                        triangulation
                    );
                    if(!status) return 0;

                    // compute local order of regions
                    status = this->computeLocalOrderOfRegions<idType>(
                        localOrders,
                        outputOrder, // used here to temporarily store distance field

                        triangulation,
                        regionMask,
                        inputOrders,
                        masterPropagations,
                        useRegionBasedIterations
                    );
                    if(!status) return 0;

                } else {

                    // sort propagations
                    status = this->sortPropagations<idType>(
                        propagations,

                        inputOrders
                    );
                    if(!status) return 0;

                    // compute propagations
                    status = this->computeInterleaving<idType>(
                        localOrders,
                        regionMask,
                        queueMask,
                        propagationMask,
                        propagations,
                        outputOrder, // used here to temporarily store distance field

                        unauthorizedMaxima,
                        triangulation,
                        inputOrders,
                        useRegionBasedIterations
                    );
                    if(!status) return 0;

                    // finalize master propagations
                    status = this->finalizePropagations<idType>(
                        masterPropagations,
                        propagations,

                        nVertices
                    );
                    if(!status) return 0;
                    nRemovedMaxima = masterPropagations.size();
                }

                // flatten orders
                status = this->flattenOrders<idType>(
                    outputOrder,

                    masterPropagations,
                    inputOrders,
                    nVertices
                );
                if(!status) return 0;

                // compute global orders
                status = this->computeGlobalOrder<idType,idType>(
                    outputOrder,
                    sortedIndices,

                    outputOrder,
                    localOrders,
                    nVertices
                );
                if(!status) return 0;

                return 1;
            };

            template<typename idType, typename dataType>
            int detectAndRemoveMaximaByPersistence(
                std::vector<idType>& maxima,
                idType* outputOrder,
                idType* localOrders,
                idType* regionMask,
                idType* queueMask,
                Propagation<idType>** propagationMask,
                std::vector<Propagation<idType>>& propagations,
                std::vector<Propagation<idType>*>& masterPropagations,
                std::vector<std::tuple<idType,idType,idType>>& sortedIndices,
                idType& nRemovedMaxima,

                const ttk::Triangulation* triangulation,
                const dataType& persistenceThreshold,
                const dataType* inputScalars,
                const idType* inputOrders,
                const bool& useRegionBasedIterations,
                const idType& escapeInterval
            ) const {
                const idType nVertices = triangulation->getNumberOfVertices();

                int status = 0;

                // Classify Critical Points
                status = this->detectMaxima<idType>(
                    maxima,

                    triangulation,
                    inputOrders
                );
                if(!status) return 0;

                // sort critical points
                status = this->sortMaxima<idType>(
                    maxima,

                    inputOrders
                );
                if(!status) return 0;

                // init propagations
                status = this->initializePropagations<idType>(
                    propagations,
                    regionMask,
                    queueMask,
                    propagationMask,

                    maxima,
                    nVertices,

                    localOrders
                );
                if(!status) return 0;

                // compute propagations
                status = this->computePersistentPropagations<idType,dataType>(
                    regionMask, // used here as saddle mask
                    queueMask,
                    propagationMask,
                    propagations,

                    triangulation,
                    inputOrders,
                    inputScalars,
                    persistenceThreshold,
                    escapeInterval
                );

                // compute trunk
                status = this->computeTrunkII<idType>(
                    propagationMask,
                    propagations,

                    triangulation,
                    regionMask, // used here as saddle mask
                    inputOrders,
                    sortedIndices,
                    inputScalars,
                    persistenceThreshold
                );
                if(!status) return 0;

                // finalize master propagations
                status = this->finalizePropagationsByPersistence<idType,dataType>(
                    masterPropagations,
                    propagations,

                    inputScalars,
                    nVertices
                );
                if(!status) return 0;

                // compute regions
                status = this->computeRegions<idType>(
                    regionMask,
                    propagationMask,
                    masterPropagations,

                    triangulation
                );
                if(!status) return 0;

                // compute local order of regions
                status = this->computeLocalOrderOfRegions<idType>(
                    localOrders,
                    outputOrder, // used here to temporarily store distance field

                    triangulation,
                    regionMask,
                    inputOrders,
                    masterPropagations,
                    useRegionBasedIterations
                );
                if(!status) return 0;

                // flatten orders
                status = this->flattenOrders<idType>(
                    outputOrder,

                    masterPropagations,
                    inputOrders,
                    nVertices
                );
                if(!status) return 0;

                // compute global orders
                status = this->computeGlobalOrder<idType,idType>(
                    outputOrder,
                    sortedIndices,

                    outputOrder,
                    localOrders,
                    nVertices
                );
                if(!status) return 0;

                nRemovedMaxima = masterPropagations.size();

                return 1;
            };

            template<typename idType, typename dataType>
            int allocateMemory(
                std::vector<idType>* inputOrders,
                std::vector<idType>* unauthorizedExtrema,
                std::vector<idType>* regionMask,
                std::vector<idType>* queueMask,
                std::vector<Propagation<idType>*>* propagationMask,
                std::vector<idType>* localOrders,
                std::vector<std::tuple<idType,idType,idType>>* sortedIndices,
                std::vector<std::tuple<dataType,idType,idType>>* sortedIndicesII,
                std::vector<Propagation<idType>>* propagationsMax,
                std::vector<Propagation<idType>*>* masterPropagationsMax,
                std::vector<Propagation<idType>>* propagationsMin,
                std::vector<Propagation<idType>*>* masterPropagationsMin,

                const idType& nVertices
            ) const {
                ttk::Timer timer;

                // allocate memory
                this->printMsg(
                    "Allocating memory",
                    0,0,this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                if(inputOrders){
                    inputOrders->resize(nVertices);
                    auto& inputOrders_ = *inputOrders;
                    #pragma omp parallel for num_threads(this->threadNumber_)
                    for(idType v=0; v<nVertices; v++)
                        inputOrders_[v] = v;
                }
                if(unauthorizedExtrema) unauthorizedExtrema->resize(nVertices);
                if(regionMask) regionMask->resize(nVertices);
                if(queueMask) queueMask->resize(nVertices);
                if(propagationMask) propagationMask->resize(nVertices);
                if(localOrders) localOrders->resize(nVertices);
                if(sortedIndices) sortedIndices->resize(nVertices);
                if(sortedIndicesII) sortedIndicesII->resize(nVertices);

                if(propagationsMax) propagationsMax->clear();
                if(masterPropagationsMax) masterPropagationsMax->clear();
                if(propagationsMin) propagationsMin->clear();
                if(masterPropagationsMin) masterPropagationsMin->clear();

                this->printMsg(
                    "Allocating memory",
                    1,timer.getElapsedTime(),this->threadNumber_
                );
                this->printMsg(debug::Separator::L2);

                return 1;
            }

            template<typename t>
            int deallocateMemory_(
                std::vector<t>* vector
            ) const {
                if(vector==nullptr)
                    return 0;

                #pragma omp task
                {
                    vector->clear();
                    auto temp = std::vector<t>();
                    (*vector) = temp;
                }
                return 1;
            }


            template<typename idType, typename dataType>
            int deallocateMemory(
                std::vector<idType>* inputOrders,
                std::vector<idType>* unauthorizedExtrema,
                std::vector<idType>* regionMask,
                std::vector<idType>* queueMask,
                std::vector<Propagation<idType>*>* propagationMask,
                std::vector<idType>* localOrders,
                std::vector<std::tuple<idType,idType,idType>>* sortedIndices,
                std::vector<std::tuple<dataType,idType,idType>>* sortedIndicesII,
                std::vector<Propagation<idType>>* propagationsMax,
                std::vector<Propagation<idType>*>* masterPropagationsMax,
                std::vector<Propagation<idType>>* propagationsMin,
                std::vector<Propagation<idType>*>* masterPropagationsMin
            ) const {
                ttk::Timer timer;

                this->printMsg(debug::Separator::L2);
                this->printMsg(
                    "Deallocating memory",
                    0,0,this->threadNumber_,
                    debug::LineMode::REPLACE
                );


                #pragma omp parallel num_threads(this->threadNumber_)
                #pragma omp single
                {
                    this->deallocateMemory_(propagationsMin);
                    this->deallocateMemory_(propagationsMax);
                    this->deallocateMemory_(inputOrders);
                    this->deallocateMemory_(unauthorizedExtrema);
                    this->deallocateMemory_(regionMask);
                    this->deallocateMemory_(queueMask);
                    this->deallocateMemory_(propagationMask);
                    this->deallocateMemory_(localOrders);
                    this->deallocateMemory_(sortedIndices);
                    this->deallocateMemory_(sortedIndicesII);
                    this->deallocateMemory_(masterPropagationsMax);
                    this->deallocateMemory_(masterPropagationsMin);
                }

                this->printMsg(
                    "Deallocating memory",
                    1,timer.getElapsedTime(),this->threadNumber_
                );

                return 1;
            }

            template<typename dataType, typename idType>
            int removeUnauthorizedExtrema(
                dataType* outputScalars,
                idType* outputOrder,

                const ttk::Triangulation* triangulation,
                const dataType* inputScalars,
                const idType* authorizedExtremaIndices,
                const idType& nAuthorizedExtremaIndices,
                const bool&   useRegionBasedIterations,
                const bool&   useInterleaving,
                const bool&   enforceAuthorizedExtrema,
                const bool&   addPerturbation,
                const bool&   useExplicitDeallocation
            ) const {

                this->printMsg(debug::Separator::L1);
                this->printMsg({
                    {"Use Region-Based Iterations", std::string(useRegionBasedIterations ? "true" : "false")},
                    {"Use Interleaving", std::string(useInterleaving ? "true" : "false")},
                    {"Enforce Authorized Extrema", std::string(enforceAuthorizedExtrema ? "true" : "false")},
                    {"Add Perturbation", std::string(addPerturbation ? "true" : "false")},
                    {"Use Explicit Deallocation", std::string(useExplicitDeallocation ? "true" : "false")}
                });
                this->printMsg(debug::Separator::L2);

                idType nVertices = triangulation->getNumberOfVertices();

                // Allocating Memory
                std::vector<idType> inputOrders;
                std::vector<idType> unauthorizedExtrema;
                std::vector<idType> regionMask;
                std::vector<idType> queueMask;
                std::vector<Propagation<idType>*> propagationMask;
                std::vector<idType> localOrders;
                std::vector<std::tuple<idType,idType,idType>> sortedIndices;
                std::vector<std::tuple<dataType,idType,idType>> sortedIndicesII;
                std::vector<Propagation<idType>> propagationsMax;
                std::vector<Propagation<idType>*> masterPropagationsMax;
                std::vector<Propagation<idType>> propagationsMin;
                std::vector<Propagation<idType>*> masterPropagationsMin;

                this->allocateMemory<idType>(
                    &inputOrders,
                    &unauthorizedExtrema,
                    &regionMask,
                    &queueMask,
                    &propagationMask,
                    &localOrders,
                    &sortedIndices,
                    &sortedIndicesII,
                    &propagationsMax,
                    &masterPropagationsMax,
                    &propagationsMin,
                    &masterPropagationsMin,

                    nVertices
                );

                // Initialize orders and scalars
                ttk::Timer timer;
                int status = 0;
                {
                    status = this->computeGlobalOrder<dataType,idType>(
                        outputOrder,
                        sortedIndicesII,

                        inputScalars,
                        inputOrders.data(),
                        nVertices
                    );
                    if(!status) return 0;

                    status = this->initializeScalars<idType,dataType>(
                        outputScalars,

                        inputScalars,
                        nVertices
                    );
                    if(!status) return 0;
                }

                // execute iterations
                size_t iteration=0;
                int sortDirection = 0;
                while(true){
                    if(!useRegionBasedIterations)
                        this->printMsg(
                            "Iteration: "+std::to_string(iteration++),
                            ttk::debug::Separator::L2
                        );

                    idType nRemovedMinima=0;
                    idType nRemovedMaxima=0;

                    // Minima
                    {
                        this->printMsg("----------- [Removing unauthorized minima]", ttk::debug::Separator::L2);

                        // invert orders to first remove minima (now maxima)
                        status = this->invertField<idType>(
                            outputOrder,
                            inputOrders.data(),

                            nVertices
                        );
                        if(!status) return 0;

                        status = this->detectAndRemoveUnauthorizedMaxima<idType>(
                            unauthorizedExtrema,
                            outputOrder,
                            localOrders.data(),
                            regionMask.data(),
                            queueMask.data(),
                            propagationMask.data(),
                            propagationsMin,
                            masterPropagationsMin,
                            sortedIndices,
                            nRemovedMinima,

                            triangulation,
                            authorizedExtremaIndices,
                            nAuthorizedExtremaIndices,
                            inputOrders.data(),
                            useRegionBasedIterations,
                            useInterleaving
                        );
                        if(!status) return 0;

                        if(nRemovedMinima){
                            sortDirection=-1;
                            status = this->flattenScalars<dataType,idType>(
                                outputScalars,

                                masterPropagationsMin
                            );
                            if(!status) return 0;
                        }
                    }

                    // Maxima
                    {
                        this->printMsg("----------- [Removing unauthorized maxima]", ttk::debug::Separator::L2);

                        // invert orders again to now remove maxima
                        status = this->invertField<idType>(
                            outputOrder,
                            inputOrders.data(),

                            nVertices
                        );
                        if(!status) return 0;

                        status = this->detectAndRemoveUnauthorizedMaxima<idType>(
                            unauthorizedExtrema,
                            outputOrder,
                            localOrders.data(),
                            regionMask.data(),
                            queueMask.data(),
                            propagationMask.data(),
                            propagationsMax,
                            masterPropagationsMax,
                            sortedIndices,
                            nRemovedMaxima,

                            triangulation,
                            authorizedExtremaIndices,
                            nAuthorizedExtremaIndices,
                            inputOrders.data(),
                            useRegionBasedIterations,
                            useInterleaving
                        );
                        if(!status) return 0;

                        if(nRemovedMaxima){
                            sortDirection=+1;
                            status = this->flattenScalars<dataType,idType>(
                                outputScalars,

                                masterPropagationsMax
                            );
                            if(!status) return 0;
                        }
                    }

                    if(useRegionBasedIterations || (nRemovedMinima+nRemovedMaxima)==0)
                        break;
                }

                if(enforceAuthorizedExtrema){
                    status = this->enforceAuthorizedExtrema<idType,dataType>(
                        outputOrder,

                        triangulation,
                        inputScalars,
                        authorizedExtremaIndices,
                        nAuthorizedExtremaIndices
                    );
                    if(!status) return 0;
                }

                // optionally add perturbation
                if(addPerturbation && sortDirection!=0){
                    this->printMsg(debug::Separator::L2);
                    status = this->computeNumericalPerturbation<dataType,idType>(
                        outputScalars,

                        outputOrder,
                        sortedIndices,
                        sortDirection
                    );
                    if(!status) return 0;
                }

                this->printMsg(debug::Separator::L2);
                this->printMsg("Complete", 1, timer.getElapsedTime(), this->threadNumber_);

                if(useExplicitDeallocation){
                    status = this->deallocateMemory<idType,dataType>(
                        &inputOrders,
                        &unauthorizedExtrema,
                        &regionMask,
                        &queueMask,
                        &propagationMask,
                        &localOrders,
                        &sortedIndices,
                        &sortedIndicesII,
                        &propagationsMax,
                        &masterPropagationsMax,
                        &propagationsMin,
                        &masterPropagationsMin
                    );
                    if(!status) return 0;
                }

                this->printMsg(debug::Separator::L1);

                return 1;
            };

            template<typename idType, typename dataType>
            int removeExtremaByPersistence(
                dataType* outputScalars,
                idType* outputOrder,

                const ttk::Triangulation* triangulation,
                const dataType* inputScalars,
                const dataType& persistenceThreshold,
                const bool& useRegionBasedIterations,
                const bool& addPerturbation,
                const bool& useExplicitDeallocation,
                const idType& escapeInterval = 1000
            ) const {

                this->printMsg(debug::Separator::L1);
                this->printMsg({
                    {"Persistence Threshold", std::to_string(persistenceThreshold)},
                    {"Escape Interval", std::to_string(escapeInterval)},
                    {"Use Region-Based Iterations", std::string(useRegionBasedIterations ? "true" : "false")},
                    {"Add Perturbation", std::string(addPerturbation ? "true" : "false")},
                    {"Use Explicit Deallocation", std::string(useExplicitDeallocation ? "true" : "false")}
                });
                this->printMsg(debug::Separator::L2);

                idType nVertices = triangulation->getNumberOfVertices();

                // Allocating Memory
                std::vector<idType> inputOrders;
                std::vector<idType> extrema;
                std::vector<idType> regionMask;
                std::vector<idType> queueMask;
                std::vector<Propagation<idType>*> propagationMask;
                std::vector<idType> localOrders;
                std::vector<std::tuple<idType,idType,idType>> sortedIndices;
                std::vector<std::tuple<dataType,idType,idType>> sortedIndicesII;
                std::vector<Propagation<idType>> propagationsMax;
                std::vector<Propagation<idType>*> masterPropagationsMax;
                std::vector<Propagation<idType>> propagationsMin;
                std::vector<Propagation<idType>*> masterPropagationsMin;

                this->allocateMemory<idType>(
                    &inputOrders,
                    &extrema,
                    &regionMask,
                    &queueMask,
                    &propagationMask,
                    &localOrders,
                    &sortedIndices,
                    &sortedIndicesII,
                    &propagationsMax,
                    &masterPropagationsMax,
                    &propagationsMin,
                    &masterPropagationsMin,

                    nVertices
                );

                // Initialize orders and scalars
                ttk::Timer timer;
                int status = 0;
                {
                    status = this->computeGlobalOrder<dataType,idType>(
                        outputOrder,
                        sortedIndicesII,

                        inputScalars,
                        inputOrders.data(),
                        nVertices,

                        sortedIndices.data()
                    );
                    if(!status) return 0;

                    status = this->initializeScalars<idType,dataType>(
                        outputScalars,

                        inputScalars,
                        nVertices
                    );
                    if(!status) return 0;
                }

                // execute iterations
                size_t iteration=0;
                int sortDirection = 0;
                while(true){

                    if(!useRegionBasedIterations)
                        this->printMsg(
                            "Iteration: "+std::to_string(iteration++),
                            ttk::debug::Separator::L2
                        );

                    idType nRemovedMinima=0;
                    idType nRemovedMaxima=0;

                    // Minima
                    {
                        this->printMsg("----------- [Removing non-persistent minima]", ttk::debug::Separator::L2);

                        // invert orders to first remove minima (now maxima)
                        status = this->invertField<idType>(
                            outputOrder,
                            inputOrders.data(),

                            nVertices
                        );
                        if(!status) return 0;

                        status = this->detectAndRemoveMaximaByPersistence<idType,dataType>(
                            extrema,
                            outputOrder,
                            localOrders.data(),
                            regionMask.data(),
                            queueMask.data(),
                            propagationMask.data(),
                            propagationsMin,
                            masterPropagationsMin,
                            sortedIndices,
                            nRemovedMinima,

                            triangulation,
                            persistenceThreshold,
                            outputScalars,
                            inputOrders.data(),
                            useRegionBasedIterations,
                            escapeInterval
                        );
                        if(!status) return 0;

                        if(nRemovedMinima){
                            sortDirection=-1;
                            status = this->flattenScalars<dataType,idType>(
                                outputScalars,

                                masterPropagationsMin
                            );
                            if(!status) return 0;
                        }
                    }

                    // Maxima
                    {
                        this->printMsg("----------- [Removing non-persistent maxima]", ttk::debug::Separator::L2);

                        // invert orders again to now remove maxima
                        status = this->invertField<idType>(
                            outputOrder,
                            inputOrders.data(),

                            nVertices
                        );
                        if(!status) return 0;

                        status = this->detectAndRemoveMaximaByPersistence<idType,dataType>(
                            extrema,
                            outputOrder,
                            localOrders.data(),
                            regionMask.data(),
                            queueMask.data(),
                            propagationMask.data(),
                            propagationsMax,
                            masterPropagationsMax,
                            sortedIndices,
                            nRemovedMaxima,

                            triangulation,
                            persistenceThreshold,
                            outputScalars,
                            inputOrders.data(),
                            useRegionBasedIterations,
                            escapeInterval
                        );
                        if(!status) return 0;

                        if(nRemovedMaxima){
                            sortDirection=+1;
                            status = this->flattenScalars<dataType,idType>(
                                outputScalars,

                                masterPropagationsMax
                            );
                            if(!status) return 0;
                        }
                    }

                    if(useRegionBasedIterations || (nRemovedMinima+nRemovedMaxima)==0)
                        break;
                }

                // optionally add perturbation
                if(addPerturbation && sortDirection!=0){
                    this->printMsg(debug::Separator::L2);
                    this->computeNumericalPerturbation<dataType,idType>(
                        outputScalars,

                        outputOrder,
                        sortedIndices,
                        sortDirection
                    );
                }

                this->printMsg(debug::Separator::L2);
                this->printMsg("Complete", 1, timer.getElapsedTime(), this->threadNumber_);

                if(useExplicitDeallocation){
                    this->deallocateMemory<idType,dataType>(
                        &inputOrders,
                        &extrema,
                        &regionMask,
                        &queueMask,
                        &propagationMask,
                        &localOrders,
                        &sortedIndices,
                        &sortedIndicesII,
                        &propagationsMax,
                        &masterPropagationsMax,
                        &propagationsMin,
                        &masterPropagationsMin
                    );
                }

                this->printMsg(debug::Separator::L1);

                return 1;
            }
    };
}