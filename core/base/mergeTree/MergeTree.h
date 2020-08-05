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

namespace ttk {

    class MergeTree : virtual public Debug {

        public:

            MergeTree(){
                this->setDebugMsgPrefix("MergeTree");
            };
            ~MergeTree(){};

            /// TODO
            /// TODO
            /// TODO
            int PreconditionTriangulation(
                ttk::Triangulation* triangulation
            ) const {
                triangulation->preconditionVertexNeighbors();
                return 1;
            }

            /// TODO
            /// TODO
            /// TODO
            template<typename dataType, typename IT>
            int computeOrderField(
                IT* orderField,

                const IT& nVertices,
                const dataType* rank1,
                const IT* rank2 = nullptr
            ) const {
                ttk::Timer timer;

                this->printMsg(
                    "Computing Order Field",
                    0, timer.getElapsedTime(), this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                // init tuples
                std::vector<std::tuple<dataType,IT,IT>> indices(nVertices);
                if(rank2==nullptr){
                    #pragma omp parallel for num_threads(this->threadNumber_)
                    for(IT i=0; i<nVertices; i++){
                        auto& t = indices[i];
                        std::get<0>(t) = rank1[i];
                        std::get<1>(t) = i;
                        std::get<2>(t) = i;
                    }
                } else {
                    #pragma omp parallel for num_threads(this->threadNumber_)
                    for(IT i=0; i<nVertices; i++){
                        auto& t = indices[i];
                        std::get<0>(t) = rank1[i];
                        std::get<1>(t) = rank2[i];
                        std::get<2>(t) = i;
                    }
                }

                this->printMsg(
                    "Computing Order Field",
                    0.2, timer.getElapsedTime(), this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                #if TTK_ENABLE_OPENMP && !defined __clang__
                    omp_set_num_threads(this->threadNumber_);
                    __gnu_parallel::sort(indices.begin(), indices.end());
                    omp_set_num_threads(1);
                #else
                    this->printWrn("Caution, outside GCC, sequential sort");
                    std::sort(indices.begin(), indices.end());
                #endif

                this->printMsg(
                    "Computing Order Field",
                    0.8, timer.getElapsedTime(), this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                #pragma omp parallel for num_threads(this->threadNumber_)
                for(IT i=0; i<nVertices; i++)
                    orderField[std::get<2>(indices[i])] = i;

                this->printMsg(
                    "Computing Order Field",
                    1, timer.getElapsedTime(), this->threadNumber_
                );

                return 1;
            }

            /// TODO
            /// TODO
            /// TODO
            template<typename IT, typename TT>
            int initializePropagations(
                std::vector<Propagation<IT>>& propagations,

                const TT* triangulation,
                const IT* orderField
            ) const {

                ttk::Timer timer;
                this->printMsg(
                    "Initialize Propagations",
                    0,0,this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                const IT nVertices = triangulation->getNumberOfVertices();
                IT writeIndex=0;

                // find discareded maxima
                #pragma omp parallel for num_threads(this->threadNumber_)
                for(IT v=0; v<nVertices; v++){
                    // check if v has larger neighbors
                    bool hasLargerNeighbor = false;
                    const IT& vOrder = orderField[v];
                    const IT nNeighbors = triangulation->getVertexNeighborNumber( v );
                    for(IT n=0; n<nNeighbors; n++){
                        IT u;
                        triangulation->getVertexNeighbor(v,n,u);
                        if( vOrder<orderField[u] ){
                            hasLargerNeighbor = true;
                            break;
                        }
                    }

                    // if v has larger neighbors then v can not be maximum
                    if(hasLargerNeighbor)
                        continue;

                    // get local write index for this thread
                    IT localWriteIndex = 0;
                    #pragma omp atomic capture
                    localWriteIndex = writeIndex++;

                    // write maximum index
                    propagations[localWriteIndex].criticalPoints.push_back(v);
                }

                // resize propagations to correct size
                propagations.resize(writeIndex);

                this->printMsg(
                    "Initialize Propagations ("+std::to_string(writeIndex)+"|"+std::to_string(nVertices)+")",
                    1,timer.getElapsedTime(),this->threadNumber_
                );

                return 1;
            }

            /// TODO
            /// TODO
            /// TODO
            template<typename IT, typename TT>
            int computeDynamicPropagation(
                Propagation<IT>& propagation,
                IT* segmentationIds,
                Propagation<IT>** propagationMask,
                IT* queueMask,

                const TT* triangulation,
                const IT* orderField,
                const IT& nActivePropagations
            ) const {

                // pointer used to compare against representative
                auto* currentPropagation = &propagation;

                // frequently used propagation members
                IT& extremumIndex = currentPropagation->criticalPoints[0];
                auto* queue = &currentPropagation->queue;

                // add extremumIndex to queue
                queue->emplace(orderField[extremumIndex],extremumIndex);

                IT queueMaskLabel = extremumIndex;
                queueMask[extremumIndex] = queueMaskLabel;
                IT segmentationId = extremumIndex;

                IT counter = 0;

                // grow region until prop reaches a saddle and then decide if prop should continue
                IT v = -1;
                while(!queue->empty()){
                    v = std::get<1>(queue->top());
                    queue->pop();

                    // continue if this thread has already seen this vertex
                    if(propagationMask[v]!=nullptr)
                        continue;

                    const IT& orderV = orderField[v];

                    // add neighbors to queue AND check if v is a saddle
                    bool isSaddle = false;
                    const IT nNeighbors = triangulation->getVertexNeighborNumber( v );

                    IT numberOfLargerNeighbors = 0;
                    IT numberOfLargerNeighborsThisThreadVisited = 0;
                    for(IT n=0; n<nNeighbors; n++){
                        IT u;
                        triangulation->getVertexNeighbor(v,n,u);

                        const IT& orderU = orderField[u];

                        // if larger neighbor
                        if( orderU>orderV ){
                            numberOfLargerNeighbors++;

                            auto uPropagation = propagationMask[u];
                            if(uPropagation==nullptr || currentPropagation!=uPropagation->find())
                                isSaddle = true;
                            else
                                numberOfLargerNeighborsThisThreadVisited++;
                        }
                        // else {
                        //     queue->emplace(orderU,u);
                        // }
                        else if(queueMask[u] != queueMaskLabel) {
                            queue->emplace(orderU,u);
                            queueMask[u] = queueMaskLabel;
                        }
                    }

                    // if v is a saddle we have to check if the current thread is the last visitor
                    if(isSaddle){

                        currentPropagation->criticalPoints.push_back(v);

                        IT numberOfRegisteredLargerVertices=0;
                        #pragma omp atomic capture
                        {
                            segmentationIds[v] -= numberOfLargerNeighborsThisThreadVisited;
                            numberOfRegisteredLargerVertices = segmentationIds[v];
                        }

                        // if this thread did not register the last remaining larger vertices then terminate propagation
                        if(numberOfRegisteredLargerVertices != -numberOfLargerNeighbors-1)
                            return 1;

                        // merge propagations
                        for(IT n=0; n<nNeighbors; n++){
                            IT u;
                            triangulation->getVertexNeighbor(v,n,u);

                            auto uPropagation = propagationMask[u];
                            if(uPropagation!=nullptr && uPropagation->find()!=currentPropagation){
                                currentPropagation = Propagation<IT>::unify(
                                    currentPropagation,
                                    uPropagation,
                                    orderField
                                );
                            }
                        }

                        queue = &currentPropagation->queue;
                        segmentationId = v;
                        queueMaskLabel = currentPropagation->criticalPoints[0];
                    }

                    // mark vertex as visited and continue
                    propagationMask[v] = currentPropagation;
                    segmentationIds[v] = segmentationId;

                    if(counter++>100){
                        counter = 0;

                        IT nActivePropagations_;

                        #pragma omp atomic read
                        nActivePropagations_ = nActivePropagations;

                        if(nActivePropagations_==1)
                            return 1;
                    }
                }

                currentPropagation->criticalPoints.push_back(v);

                return 1;
            }

            /// TODO
            /// TODO
            /// TODO
            template<typename IT, typename TT>
            int computeDynamicPropagations(
                IT* segmentationIds,
                std::vector<Propagation<IT>>&  propagations,
                Propagation<IT>** propagationMask,
                IT* queueMask,

                const TT* triangulation,
                const IT* orderField
            ) const {
                ttk::Timer timer;

                int status = 1;
                const IT nPropagations = propagations.size();
                IT nActivePropagations = nPropagations;

                this->printMsg(
                    "Computing Dynamic Propagations ("+std::to_string(nPropagations)+")",
                    0, 0, this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                // compute propagations
                #pragma omp parallel for schedule(dynamic,1) num_threads(this->threadNumber_)
                for(IT p=0; p<nPropagations; p++){
                    int localStatus = this->computeDynamicPropagation<IT,TT>(
                        propagations[p],
                        segmentationIds,
                        propagationMask,
                        queueMask,

                        triangulation,
                        orderField,
                        nActivePropagations
                    );
                    if(!localStatus)
                        status = 0;

                    #pragma omp atomic update
                    nActivePropagations--;
                }
                if(!status) return 0;

                this->printMsg(
                    "Computing Dynamic Propagations ("+std::to_string(nPropagations)+")",
                    1, timer.getElapsedTime(), this->threadNumber_
                );

                return 1;
            }

            template<typename IT, typename TT>
            int computeTrunk(
                IT* segmentationIds,
                std::vector<Propagation<IT>>& propagations,
                Propagation<IT>** propagationMask,
                IT* mask,

                const TT* triangulation,
                const IT* orderField
            ) const {
                ttk::Timer timer;
                this->printMsg(
                    "Computing trunk",
                    0, 0, this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                const IT nVertices = triangulation->getNumberOfVertices();

                // reconsstruct sorted vertices
                for(IT v=0; v<nVertices; v++)
                    mask[orderField[v]] = v;

                // find largest propagation
                Propagation<IT>* currentPropagation = nullptr;
                IT segmentationId = -1;
                IT trunkIndex = nVertices-1;
                for(; trunkIndex>=0; trunkIndex--){
                    const IT& v = mask[trunkIndex];

                    if(propagationMask[v]!=nullptr)
                        continue;

                    IT nNeighbors = triangulation->getVertexNeighborNumber(v);
                    for(IT n=0; n<nNeighbors; n++){
                        IT u;
                        triangulation->getVertexNeighbor(v,n,u);

                        if(propagationMask[u]!=nullptr){
                            currentPropagation = propagationMask[u]->find();
                            segmentationId = segmentationIds[u];
                            break;
                        }
                    }

                    break;
                }

                if(currentPropagation==nullptr){
                    this->printWrn("Empty Trunk.");
                    return 1;
                }

                IT nTrunkVertices = 0;
                for(; trunkIndex>=0; trunkIndex--){
                    const IT& v = mask[trunkIndex];

                    if(propagationMask[v]!=nullptr)
                        continue;

                    nTrunkVertices++;

                    // if saddle
                    if(segmentationIds[v]<-1){

                        currentPropagation->criticalPoints.push_back(v);
                        const IT nNeighbors = triangulation->getVertexNeighborNumber( v );

                        // merge propagations
                        for(IT n=0; n<nNeighbors; n++){
                            IT u;
                            triangulation->getVertexNeighbor(v,n,u);

                            auto uPropagation = propagationMask[u];
                            if(uPropagation!=nullptr && uPropagation->find()!=currentPropagation){
                                currentPropagation = Propagation<IT>::unify(
                                    currentPropagation,
                                    uPropagation,
                                    orderField
                                );
                            }
                        }

                        segmentationId = v;
                    }

                    propagationMask[v] = currentPropagation;
                    segmentationIds[v] = segmentationId;
                }

                // add global minimum
                currentPropagation->criticalPoints.push_back( mask[0] );

                std::stringstream vFraction;
                vFraction << std::fixed << std::setprecision(2) << ((float)nTrunkVertices/(float)nVertices);

                this->printMsg(
                    "Computing trunk ("+std::to_string(nTrunkVertices)+"|"+vFraction.str()+")",
                    1, timer.getElapsedTime(), this->threadNumber_
                );

                return 1;
            }

            /// TODO
            /// TODO
            /// TODO
            template<typename IT>
            int allocateMemory(
                std::vector<Propagation<IT>>& propagations,
                std::vector<Propagation<IT>*>& propagationMask,
                std::vector<IT>& queueMask,
                IT* segmentationIds,

                const IT& nVertices
            ) const {
                ttk::Timer timer;

                // =============================================================
                // allocate and init memory
                // =============================================================
                this->printMsg(
                    "Allocating Memory",
                    0,0,this->threadNumber_,
                    debug::LineMode::REPLACE
                );

                propagations.resize(nVertices);
                propagationMask.resize(nVertices, nullptr);
                queueMask.resize(nVertices, -1);

                #pragma omp parallel for num_threads(this->threadNumber_)
                for(IT i=0; i<nVertices; i++)
                    segmentationIds[i] = -1;

                this->printMsg(
                    "Allocating Memory",
                    1,timer.getElapsedTime(),this->threadNumber_
                );

                return 1;
            }

            /// TODO
            /// TODO
            /// TODO
            template<typename IT, typename TT>
            int computeMergeTreeSegmentation(
                IT* segmentationIds,
                std::vector<Propagation<IT>>& propagations,

                const TT* triangulation,
                const IT* orderField
            ) const {
                ttk::Timer timer;

                this->printMsg(debug::Separator::L2);

                const IT nVertices = triangulation->getNumberOfVertices();

                // =============================================================
                // allocate and init memory
                // =============================================================
                std::vector<Propagation<IT>*> propagationMask;
                std::vector<IT> queueMask;
                if(!this->allocateMemory<IT>(
                    propagations,
                    propagationMask,
                    queueMask,
                    segmentationIds,

                    nVertices
                ))
                    return 0;

                // =============================================================
                // initialize propagations
                // =============================================================
                if(!this->initializePropagations<IT, TT>(
                    propagations,

                    triangulation,
                    orderField
                ))
                    return 0;

                // =============================================================
                // execute propagations
                // =============================================================
                if(!this->computeDynamicPropagations<IT, TT>(
                    segmentationIds,
                    propagations,
                    propagationMask.data(),
                    queueMask.data(),

                    triangulation,
                    orderField
                ))
                    return 0;

                // =============================================================
                // compute trunk
                // =============================================================
                if(!this->computeTrunk<IT, TT>(
                    segmentationIds,
                    propagations,
                    propagationMask.data(),
                    queueMask.data(),

                    triangulation,
                    orderField
                ))
                    return 0;

                this->printMsg(debug::Separator::L2);
                this->printMsg(
                    "Complete",
                    1,timer.getElapsedTime(),this->threadNumber_
                );
                this->printMsg(debug::Separator::L2);

                return 1;
            }
    };
}