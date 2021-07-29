/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::TrackingGraph
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// This module defines the %TrackingGraph class that computes for each vertex
/// of a triangulation the average scalar value of itself and its direct
/// neighbors.
///
/// \b Related \b publication: \n
/// 'TrackingGraph'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

#include <unordered_map>
#include <unordered_set>

namespace ttk {

  class TrackingGraph : virtual public Debug {

  public:

    struct Edge {
      int u;
      int v;
      int e;
    };

    std::vector<std::vector<Edge>> inEdges;
    std::vector<std::vector<Edge>> outEdges;

    TrackingGraph();

    template<typename IT>
    int preconditionInOutEdges(
      const int nNodes,
      const int nEdges,
      const IT* connectivityList
    ) {
      ttk::Timer timer;
      const std::string msg("Building Tracking Graph Structure");
      this->printMsg(msg,0,0,1,ttk::debug::LineMode::REPLACE);

      this->outEdges.clear();
      this->outEdges.resize(nNodes);
      this->inEdges.clear();
      this->inEdges.resize(nNodes);

      // std::vector<std::unordered_set<int>> nextNodesTemp;
      // std::vector<std::unordered_set<int>> prevNodesTemp;
      // nextNodesTemp.resize(nNodes);
      // prevNodesTemp.resize(nNodes);

      for(int i=0; i<nEdges; i++){
        const int u = connectivityList[i*2+0];
        const int v = connectivityList[i*2+1];

        this->outEdges[u].push_back({u,v,i});
        this->inEdges[v].push_back({u,v,i});
      }

      // for(int i=0; i<nNodes; i++){
      //   this->nextNodes[i].insert(this->nextNodes[i].begin(),nextNodesTemp[i].begin(),nextNodesTemp[i].end());
      //   this->prevNodes[i].insert(this->prevNodes[i].begin(),prevNodesTemp[i].begin(),prevNodesTemp[i].end());
      // }

      this->printMsg(msg,1,timer.getElapsedTime(),1);
      return 1;
    }

    // template<typename IT>
    // int preconditionPrevAndNextNodesByTime(
    //   const int nNodes,
    //   const int nEdges,
    //   const IT* connectivityList,
    //   const int* time
    // ) {
    //   ttk::Timer timer;
    //   const std::string msg("Building Tracking Graph Structure");
    //   this->printMsg(msg,0,0,1,ttk::debug::LineMode::REPLACE);

    //   this->prevNodesByTime.clear();
    //   this->nextNodesByTime.clear();

    //   for(int i=0; i<nNodes; i++){
    //     this->prevNodesByTime.insert({time[i],std::unordered_map<int,std::unordered_set<int>>()});
    //     this->nextNodesByTime.insert({time[i],std::unordered_map<int,std::unordered_set<int>>()});
    //   }

    //   for(int i=0; i<nEdges; i++){
    //     const int u = connectivityList[i*2+0];
    //     const int v = connectivityList[i*2+1];
    //     this->nextNodesByTime.find(time[u])->second.find(u)->second.insert(v);
    //     this->prevNodesByTime.find(time[v])->second.find(v)->second.insert(u);
    //   }

    //   this->printMsg(msg,1,timer.getElapsedTime(),1);
    //   return 1;
    // }
  }; // TrackingGraph class

} // namespace ttk
