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

#include "nodeContainer.hpp"
#include <valarray>
#include <fstream>
#include "exceptions.hpp"

#ifndef GRAPH
#define GRAPH

/*! \brief Network class*/

class GraphReader{
  friend class GraphBuilder;
  friend class TestGraphReader;
    // Members
    string net_str;
    vector < string > node_labels;
    vector < string > subTreeStrs;
    vector < string > brchlens;

    vector <string> tip_name;
    vector <string> tax_name;
    vector <string> sampleNames;

    // Methods
    GraphReader( string net_str );
    ~GraphReader(){ }

    void check_Parenthesis(string &in_str);
    void check_labeled( string in_str );
    void extract_tax_and_tip_names();
    void extractSampleNames();
    size_t Parenthesis_balance_index_backwards( string &in_str, size_t i );
    string label_interior_node(string in_str);
    string extract_One_subTreeStr( string &in_str, size_t back_parenthesis_index );
    string extract_label(string &in_str, size_t i);
};


bool start_of_tax_name(string in_str, size_t i);


class GraphBuilder{
  friend class TestGraphBuild;
  friend class Net;
  friend class CoalST;
  friend class CoalSN;
  friend class CoalGT;
  friend class TmpSN;
  friend class HybridCoal;
  friend class Frequency;
  friend class Figure;
  private:
    // Members
    bool is_ultrametric; /*!< \brief true if the distances between tips and root are equal; false, otherwise */ // This is used in Figure
    bool isNet_; /*!< \brief true if Net is a network; false if it's a tree */
    NodeContainer nodes_;
    GraphReader * Tree_info;
    size_t tmpEdgeNameLabel_;
    vector <string> tip_name; // maybe don't need them actually...
    vector <string> tax_name; // maybe don't need them actually...
    vector <string> sampleNames;
    //string net_str; /*!< \brief species network string \todo this is new!!!*/
    size_t max_rank;

    // Methods
    void init();
    void enumerate_internal_branch( Node *node );
    void init_descendant(){}; // maybe don't need them actually...
    void init_node_clade(){}; // maybe don't need them actually...
    void rewrite_descendant();
    string rewrite_internal_subTreeStr( Node * node );
    void rewrite_node_clade(); // maybe don't need them actually...
    size_t Parenthesis_balance_index_forwards( string &in_str, size_t i );
    void check_isNet(); /*!< \brief To determin if a Net is network or not. \return is_Net */
    void check_isUltrametric(); /*!< \brief To determin if a Net is ultrametric or not. \return is_ultrametric */

    void removeOneChildInternalNode();
    void removeZeroChildHybridNode();
    void which_taxa_is_below();
    void which_tip_is_below();
    void which_sample_is_below();
    bool isNet() const { return this->isNet_ ; }
    void print();
    bool print_all_node_dout();
    void initialize_nodes( string &net_str );
    void remove_repeated_hybrid_node();
    void connect_graph();
    NodeContainer const *getNodes() const { return &nodes_; };

    GraphBuilder( string &net_str );
    GraphBuilder(const GraphBuilder &currentGraph);
    ~GraphBuilder(){};

  public:
    NodeContainer *nodes() { return &(this->nodes_); } /*!< \brief array of nodes */
    void rewrite_subTreeStr();
    string print_newick( Node * node );
    string reWritesubTreeStrAtRoot () {
        dout << " reWritesubTreeStrAtRoot " << endl;
        for ( auto it = nodes_.iterator(); it.good(); ++it){
            dout << (*it) << " " << (*it)->nodeName << " " << (*it)->subTreeStr <<endl;
        }
        string newTreeStr = this->rewrite_internal_subTreeStr(this->nodes_.back()) + this->nodes_.back()->nodeName + ";";
        dout << " return new tree string: " << newTreeStr << endl;
        return newTreeStr;
    }
};



string remove_interior_label(string in_str);
size_t end_of_label_or_bl( string &in_str, size_t i);
void readNextStringto( string &readto , int& argc_i, int argc_, char * const* argv_ );

string rm_one_child_root(string in_str);

#endif
