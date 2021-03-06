/*
 * hybrid-coal is used to compute gene tree probabilities given species network under coalescent process.
 *
 * Copyright (C) 2010 -- 2015 Sha (Joe) Zhu
 *
 * This file is part of hybrid-coal
 *
 * hybrid-coal is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <iomanip>  // std::setw
#include <iostream> //std::cout
#include <vector>
#include <string>
#include <cassert>
#include <valarray>

using namespace std;

//Unless compiled with options NDEBUG, we will produce a debug output using
//'dout' instead of cout and execute (expensive) assert statements.
#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif

#ifndef NODE
#define NODE
/*! \brief Node of a tree or network, it also represent the branch between this node and its parent node
 */

enum NAMETYPE { TAXA, TIP };


class Edge {
  friend class Node;
  friend class GraphBuilder;
  friend class Figure;
  friend class CoalST;
  friend class CoalGT;
  friend class CoalSN;
    // For any extenal edge (connected to a tip node), the edgeName is set as zero.
    // edgeName is used to label the interior branchs of the tree
    size_t edgeName_;
    double branchLength_;
    // Consider to use member to indicate the probability as well ...

    double bl() const { return this->branchLength_;}
    void setLength ( double bl ){ this->branchLength_ = bl; }

    size_t name() const {return this->edgeName_;}
    void setName( size_t name ) { this->edgeName_ = name; }

    void setBothNameAndLength( size_t name, double bl){
        this->setName ( name );
        this->setLength ( bl );
    }

    Edge(){
        this->edgeName_     = 0;
        this->branchLength_ = 0.0;
    }
    ~Edge(){}
};


class Node {
    friend class NodeIterator;
    friend class ConstNodeIterator;
    friend class NodeContainer;
    friend class GraphBuilder;
    friend class Net;
    friend class CoalGT;
    friend class CoalST;
    friend class CoalSN;
    friend class TmpSN;
    friend class simTree;
    friend class HybridLambda;
    friend class Figure;
  public:
    ~Node(){};
    // Getters and Setters
    //double height() const { return this->height_;} // This has no use for hybrid-coal
    //void set_height ( double h ){ this->height_ = h; }// This has no use for hybrid-coal

    Edge edge1;
    Edge edge2;
    string nodeName; /*!< \brief String label of a node, each node has unique name */
    //size_t node_index; /*!< \brief node index in the array, \todo use this more often!!!*/
    string subTreeStr; /*!< \brief node content, the subtree string at this node */
    bool isHybrid() const { return ( this->parent2() != NULL ) ;} /*!< \brief Hybrid node only, indicator of a hybrid node */

  private:
    // Members
    size_t rank_;     /*!< \brief rank of the node, tip node has rank one, the root has the highest rank */
    void setRank ( size_t rank ){ this->rank_ = rank; }
    bool visited_;

    valarray < size_t > taxa_below;
    valarray < int > tips_below;
    valarray < size_t > samples_below; // tips_below
    vector < Node* > interior_nodes_below; /*!< \brief list of pointers to its descndent interior nodes */
    vector < Node* > child; /*!< \brief list of pointers to its child nodes */
    Node* parent1_; /*!< \brief pointer to its parent node. */
    Node* previous_;
    Node* next_;

    string clade; /*!< \brief clade at this node, \todo this should be modified to a vector <string> */

    int num_descndnt; /*!< \brief number of the tip nodes, that are descendant from this node */
    //int num_descndnt_interior; /*!< \brief number of the interior nodes, that are descendant from this node \todo to be replaced by interior_nodes_below.size()? */
    size_t NumberOfInteriorNodesBelow() const { return this->interior_nodes_below.size(); }
    vector <double> path_time;
    //double height_; /*!< \brief distance to the bottom of the tree */  // This has no use for hybrid-coal

    bool isTip_; /*!< \brief Indicator of tip nodes. It's true, if it is a tip node, otherwise it is false. */
    bool isTip() const { // DEBUG, there is a bug here, the assertion fails
        //assert ( this->child.size() == 0);
        return this->isTip_ ;
    }
    void setIsTip ( bool TorF ){ this->isTip_ = TorF; }
    bool isBelowHybrid_; //bool descndnt_of_hybrid; /*!< \brief Indicator of descendant of hybrid nodes. It's true, if it is a descendant of hybrid nodes; false, otherwise. */
    bool isBelowHybrid() const { return this->isBelowHybrid_; }
    void setIsBelowHybrid( bool isBelow ){ this->isBelowHybrid_ = isBelow; }

    Node* parent2_; /*!< \brief Hybrid node only, pointer to its second parent node. */
    //double prob_to_hybrid_left; /*!< \brief Hybrid node only, the probability that a lineage goes to the left */

    string nodeClass; /*!< \brief Name of a node, this is not unique for nodes. e.g. if its label is A_1, name is A */

    //vector <size_t> Net_node_contains_gt_node1; /*!< Used while simulation, check if a Network node contains a gene tree node */
    //vector <size_t> Net_node_contains_gt_node2; /*!< Used while simulation, check if a Network node contains a gene tree node */

    //size_t edge() const {return this->edge_;}
    //void setEdge( size_t num ) { this->edge_ = num; }

    //size_t edge2() const {return this->edge2_;}
    //void setEdge2( size_t num ) { this->edge2_ = num; }

    bool visited() const { return this->visited_; }
    void set_visited ( bool TorF ){ this->visited_ = TorF; }

    size_t rank() const { return this->rank_; }

    //bool tip() const { return this->tip_bool; }

    Node* parent1() const { return this->parent1_ ; }
    void set_parent1 ( Node * node ) { this->parent1_ = node; }

    Node* parent2() const { return this->parent2_ ; }
    void set_parent2 ( Node * node ) { this->parent2_ = node; }

    // NodeIterator
    Node* previous() const { return this->previous_ ; }
    void set_previous ( Node * node ) { this->previous_ = node; }

    Node* next() const { return this->next_; }
    void set_next ( Node * node ) { this->next_ = node; }

    bool is_first() const { return( this->previous_ == NULL ); }
    bool is_last()  const { return( this->next_ == NULL ); }

    // Methods
    //Node (); /*!< \brief Initialize Node class*/
    Node ( size_t max_of_taxa,
           size_t max_of_tip,
           size_t max_of_sample,
           string nodeName,
           string content,
           double bl,
           bool tip );

    void init();
    void add_child( Node *child_node /*! pointer to the child node*/);
    void CalculateRank();
    void print( bool is_Net = false );
    bool print_dout( bool is_Net = false );
    //void find_tip();
    void findWhoIsBelowHybrid();
    //bool find_descndnt ( string &name, NAMETYPE type );

    double extract_hybrid_para(){
        size_t hash_index = this->nodeName.find('#');
        return strtod( this->nodeName.substr(hash_index+1).c_str(), NULL) ;
    }
};


#endif
