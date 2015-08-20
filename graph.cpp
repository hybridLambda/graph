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


#include "graph.hpp"

GraphReader::GraphReader ( string in_str ){
    // check & sign, this should be illigal for hybrid-Lambda,
    this->check_Parenthesis( in_str );
    this->check_labeled( in_str );

    // Try to use string iterator for this...
    size_t found_bl = this->net_str.find(':');
    for ( size_t charPosition = 1; charPosition < this->net_str.size(); ){
        if ( this->net_str[charPosition]=='e' && ( this->net_str[charPosition+1]=='-' || this->net_str[charPosition+1]=='+' ) ){
            charPosition++;
        }
        else if ( start_of_tax_name( this->net_str, charPosition) ){
            size_t str_start_index = charPosition;
            string label = extract_label( this->net_str, charPosition );
            this->node_labels.push_back(label);

            string subTreeStr = ( this->net_str[str_start_index-1]==')' ) ? extract_One_subTreeStr ( this->net_str, str_start_index-1 )
                                                                            : label ;
            subTreeStrs.push_back(subTreeStr);

            charPosition += label.size();

            string brchlen;
            if ( found_bl != string::npos ){
                if ( this->net_str[charPosition] != ':' ){
                    throw BranchLengthUnGiven (label);
                }
                assert ( this->net_str[charPosition] == ':' );
                size_t found=min(min( this->net_str.find(",",charPosition+1), this->net_str.find(")",charPosition+1)), this->net_str.size());
                brchlen = this->net_str.substr(charPosition+1,found-charPosition-1);
            }
            found_bl = this->net_str.find(":", found_bl+1);
            brchlens.push_back(brchlen);
        }
        else {
            charPosition++;
        }
    }
    this->extract_tax_and_tip_names();
}


/*! \brief Checking Parenthesis of a (extended) Newick string */
void GraphReader::check_Parenthesis( string &in_str ){
    int num_b = 0;
    for ( size_t i = 0; i < in_str.size(); i++){
        if      (in_str[i] == '(') num_b++;
        else if (in_str[i] == ')') num_b--;
        else continue;
    }
    if ( num_b != 0 ) throw ParenthesisNotBalanced(in_str);
}


string GraphReader::extract_label( string &in_str, size_t i ){
    string label( in_str.substr ( i , end_of_label_or_bl ( in_str, i ) + 1-i ) );
    return label;
}


string GraphReader::extract_One_subTreeStr( string &in_str, size_t back_parenthesis_index ){
    size_t front_parenthesis_index = Parenthesis_balance_index_backwards( in_str, back_parenthesis_index );
    return in_str.substr ( front_parenthesis_index, back_parenthesis_index - front_parenthesis_index + 1);
}


void GraphReader::check_labeled( string in_str ){
    bool labeled_bool=true;
    for ( size_t i = 0; i < in_str.size(); i++ ){
        if ( in_str[i] == ')' && i == end_of_label_or_bl(in_str, i) ){
            labeled_bool = false;
            break;
        }
    }
    this->net_str = labeled_bool ? in_str:label_interior_node(in_str);
}


/*! \brief Label interior node if the interior nodes of the tree string are not labeled */
string GraphReader::label_interior_node(string in_str /*!< input newick form string */){
    vector <string> in_str_partition;
    int interior_node_counter = 0;
    int sub_str_start_index = 0;
    size_t i = in_str.find(')');
    while ( i<in_str.size() ){
        interior_node_counter++;
        string current_string;
        size_t found_next_bracket = min(in_str.find(")", sub_str_start_index),in_str.size());
        current_string = in_str.substr(sub_str_start_index, found_next_bracket - sub_str_start_index +1);
        if ( in_str[i+1] == ';' || i == (in_str.size()-1) ){
            current_string += "root";
        }
        else {
            current_string += "Int_" + to_string( interior_node_counter );
            sub_str_start_index = i+1;
        }
        in_str_partition.push_back(current_string);

        i = in_str.find( ")", i+1 );
    }
    string out_str;
    for ( size_t i = 0; i < in_str_partition.size(); i++ )
        out_str += in_str_partition[i];
    return out_str;
}


void GraphReader::extract_tax_and_tip_names(){
    for ( size_t i = 0; i < node_labels.size(); i++ ){
        if ( node_labels[i] != subTreeStrs[i]) continue;         // skip interior nodes
        if ( node_labels[i].find("#") != string::npos ) continue;  // skip hybrid nodes
        if ( node_labels[i].find("_") > 0 ){
            string node_name = node_labels[i].substr( 0, node_labels[i].find("_"));
            bool new_tax_bool = true;
            for ( size_t tax_i = 0; tax_i < tax_name.size(); tax_i++ ){
                if ( tax_name[tax_i] == node_name ){
                    new_tax_bool = false;
                    break;
                }
            }
            if ( new_tax_bool ) tax_name.push_back( node_name );
        } else{
            tax_name.push_back( node_labels[i] );
        }
        tip_name.push_back( node_labels[i] );
    }
    sort(tax_name.begin(), tax_name.end());
    sort(tip_name.begin(), tip_name.end());
    this->extractSampleNames();
}


void GraphReader::extractSampleNames(){
    for (size_t i = 0; i < tip_name.size(); i++){
        // extract for & sign
        //vector <string> new_tip_list;    
        size_t found = min(tip_name[i].find("&"), tip_name[i].size());
        for ( size_t ii = 0; ii < tip_name[i].size(); ii++ ){
            sampleNames.push_back(tip_name[i].substr(ii, found - ii ));
            ii = found;
            found = min(tip_name[i].find("&", found+1 ), tip_name[i].size());
        }
        //for (size_t ii = 0 ; ii < new_tip_list.size(); ii++){
            //bool new_tax_bool=true;
            //string undetermin_taxname = new_tip_list[ii].substr(0,new_tip_list[ii].find("_"));
                
            //for (size_t tax_i = 0; tax_i < sampleNames.size(); tax_i++){
                //if (sampleNames[tax_i] == undetermin_taxname){
                    //new_tax_bool=false;
                    //break;
                //}
            //}
            //if (new_tax_bool){
                //sampleNames.push_back(undetermin_taxname);
            //}
        //}
    }
    sort(sampleNames.begin(), sampleNames.end());    
}


/*! \brief Construct Net object from a (extended) Newick string */
GraphBuilder::GraphBuilder(string &in_str /*! input (extended) newick form string */){
    if ( in_str.size() == 0 )
        return;

    this->init( );

    this->initialize_nodes( in_str );
    this->remove_repeated_hybrid_node();
    this->connect_graph();
    this->nodes_.back()->findWhoIsBelowHybrid();
    this->nodes_.back()->CalculateRank();
    this->max_rank = nodes_.back()->rank();
    this->enumerate_internal_branch( this->nodes_.back() );
    //this->init_descendant();
    //this->init_node_clade();
    //this->rewrite_descendant();
    this->check_isNet();
    //this->check_isUltrametric(); // The ultrametric function is not very useful here. As we dont actually need it in hybrid-coal, hybridlambda yes....
    //dout<<"Net constructed"<<endl;
}

GraphBuilder::GraphBuilder(const GraphBuilder &currentGraph){
    this->is_ultrametric = currentGraph.is_ultrametric;
    this->isNet_ = currentGraph.isNet_; /*!< \brief true if Net is a network; false if it's a tree */
    this->nodes_ = NodeContainer(*currentGraph.getNodes());
    this->Tree_info = NULL;
    size_t tmpEdgeNameLabel_;
    this->tip_name = currentGraph.tip_name; // maybe don't need them actually...
    this->tax_name = currentGraph.tax_name; // maybe don't need them actually...
    this->sampleNames = currentGraph.sampleNames;
    //string net_str; /*!< \brief species network string \todo this is new!!!*/
    this->max_rank = currentGraph.max_rank;
}


void GraphBuilder::init(){
    this->tmpEdgeNameLabel_ = 0;
    this->isNet_ = false;
    this->is_ultrametric = true;
}


void GraphBuilder::initialize_nodes( string &in_str ){
    this->Tree_info = new GraphReader ( in_str );
    for ( size_t i = 0 ; i < Tree_info->brchlens.size(); i++ ){
        bool isTip = ( Tree_info->node_labels[i] == Tree_info->subTreeStrs[i] && Tree_info->node_labels[i].find("#") == string::npos);
        Node * node = new Node ( Tree_info->tax_name.size(),
                                 Tree_info->tip_name.size(),
                                 Tree_info->sampleNames.size(),
                                 Tree_info->node_labels[i],
                                 Tree_info->subTreeStrs[i],
                                 strtod( Tree_info->brchlens[i].c_str(), NULL),
                                 isTip );
        this->nodes_.add ( node );

        if ( !isTip ) continue;

        size_t tip_i;
        for ( tip_i = 0 ; tip_i < Tree_info->tip_name.size(); tip_i++ ){
            if ( Tree_info->tip_name[tip_i] == Tree_info->node_labels[i] ) break;
        }
        this->nodes_.back()->tips_below[tip_i] = (size_t)1;

        size_t taxa_i;
        for ( taxa_i = 0 ; taxa_i < Tree_info->tax_name.size(); taxa_i++ ){
            if ( Tree_info->tax_name[taxa_i] == Tree_info->node_labels[i] ) break;
        }
        this->nodes_.back()->taxa_below[taxa_i] = (size_t)1;

        size_t sample_i;
        for ( sample_i = 0 ; sample_i < Tree_info->sampleNames.size(); sample_i++ ){
            if ( Tree_info->node_labels[i].find(Tree_info->sampleNames[sample_i]) != string::npos  ) {
                this->nodes_.back()->samples_below[sample_i] = (size_t)1;
            };
        }

    }
    this->tax_name = Tree_info->tax_name;
    this->tip_name = Tree_info->tip_name;
    this->sampleNames = Tree_info->sampleNames;
    delete Tree_info;
}


void GraphBuilder::remove_repeated_hybrid_node(){
    for ( auto it_i = nodes_.iterator(); it_i.good(); ++it_i){
        if ( (*it_i)->next() == NULL ) break;
        NodeIterator it_j ( (*it_i)->next() );
        for ( ; it_j.good(); ++it_j){
            if ( (*it_j)->next() == NULL ) break;
            if ( (*it_i)->nodeName != (*it_j)->nodeName ) continue;
            dout << "Remove " << (*it_j)->nodeName <<" "<< (*it_j)->subTreeStr <<endl;
            if ( (*it_j)->subTreeStr[0] == '(' ){
                (*it_i)->subTreeStr = (*it_j)->subTreeStr;
            }
            (*it_i)->edge2.setLength ( (*it_j)->edge1.bl() );
            break;
        }
        if ( (*it_j)->nodeName == (*it_i)->nodeName ) this->nodes_.remove( (*it_j) );
    }
}




void GraphBuilder::connect_graph(){
    for ( auto it = nodes_.iterator(); it.good(); ++it){
        //dout << " node " << (*it) << ": " << (*it)->nodeName <<"\t"<<(*it)->subTreeStr<<endl;
        if ( (*it)->subTreeStr[0] != '(' ) continue;

        char child_node1[(*it)->subTreeStr.length()];
        for ( size_t i_content_len = 1; i_content_len < (*it)->subTreeStr.length(); ){
            if ((*it)->subTreeStr[i_content_len]=='(' ||  start_of_tax_name((*it)->subTreeStr,i_content_len) ){
                size_t j_content_len = ((*it)->subTreeStr[i_content_len] == '(') ? Parenthesis_balance_index_forwards( (*it)->subTreeStr, i_content_len ) + 1:
                                                                                     i_content_len;
                int child1_subTreeStr_i = 0;
                for ( ; j_content_len < (*it)->subTreeStr.length(); j_content_len++){
                    child_node1[child1_subTreeStr_i] = (*it)->subTreeStr[j_content_len];
                    char stop = (*it)->subTreeStr[j_content_len+1];
                    if ( stop == ',' || stop == ')' || stop == ':'){
                        child_node1[child1_subTreeStr_i+1]='\0';
                        break;
                    }
                    child1_subTreeStr_i++;
                }
                string child_node1_str = child_node1;
                i_content_len = j_content_len + 2;
                for ( auto it_i = nodes_.iterator(); it_i.good(); ++it_i){
                        if (child_node1_str == (*it_i)->nodeName) {
                        (*it)->add_child( (*it_i) );
                        //dout << "node " << &this->nodes_[i] << " has child "<< &this->nodes_[j]<<endl;
                    }
                }
                //for ( size_t j = 0; j < nodes_.size(); j++){
                    //if (child_node1_str == this->nodes_.at(j)->label) {
                        //(*it)->add_child( this->nodes_.at(j) );
                        ////dout << "node " << &this->nodes_[i] << " has child "<< &this->nodes_[j]<<endl;
                    //}
                //}
            }
            else { i_content_len++;}
        }
    }
}






void GraphBuilder::check_isNet(){ //false stands for tree, true stands for net_work
    for ( auto it = nodes_.iterator(); it.good(); ++it){
        if ( (*it)->parent2() == NULL ) continue;
        this->isNet_ = true;
        return;
    }
}


void GraphBuilder::print(){
    if ( this->isNet() ) cout<<"           label  hybrid hyb_des non-tp parent1  abs_t brchln1 parent2 brchln2 #child #dsndnt #id rank   edge   Clade "<<endl;
    else cout<<"            label non-tp   parent       brchln #child #dsndnt #id rank edge   Clade "<<endl;
    for ( auto it = nodes_.iterator(); it.good(); ++it){
        //for (size_t j = 0; j < this->descndnt[i].size(); j++ ) {cout<<setw(3)<<this->descndnt[i][j];}
        (*it)->print( this->isNet() );
        cout<<"  ";
        //for (size_t j=0;j<this->samples_below[i].size();j++) {cout<<this->samples_below[i][j]; }
        cout<<endl;
    }
}


/*! \brief enumerate the internal branches */
void GraphBuilder::enumerate_internal_branch( Node * node ) {
    if ( node->isTip() ) return;

    if ( node->visited() ){
        this->tmpEdgeNameLabel_ ++;
        node->edge2.setName( tmpEdgeNameLabel_ );
    } else{
        for ( size_t i = 0; i < node->child.size(); i++ ){
            this->enumerate_internal_branch( node->child[i] );
        }
        node->set_visited( true );
        this->tmpEdgeNameLabel_ ++;
        node->edge1.setName( tmpEdgeNameLabel_ );
    }
}




size_t GraphReader::Parenthesis_balance_index_backwards( string &in_str, size_t i ){
    size_t j = i;
    int num_b = 0;
    for ( ; j > 0 ; j-- ){
        if      ( in_str[j] == '(' ) num_b--;
        else if ( in_str[j] == ')' ) num_b++;
        else continue;
        if ( num_b == 0 ) break;
    }
    return j;
}


size_t GraphBuilder::Parenthesis_balance_index_forwards( string &in_str, size_t i ){
    size_t j = i;
    int num_b = 0;
    for ( ; j < in_str.size(); j++ ){
        if      ( in_str[j] == '(' ) num_b++;
        else if ( in_str[j] == ')' ) num_b--;
        else continue;
        if ( num_b == 0 ) break;
    }
    return j;
}


void GraphBuilder::which_taxa_is_below(){
    for ( auto it = nodes_.iterator(); it.good(); ++it){
        Node* current_node = (*it);
        if ( !current_node->isTip() ) continue;

        size_t taxa_i = 0;
        for ( taxa_i = 0 ; taxa_i < current_node->taxa_below.size(); taxa_i++){
            if ( current_node->taxa_below[taxa_i] == 1 ) break;
        }

        while ( current_node->parent1() ){
            current_node->parent1()->taxa_below[taxa_i] = 1;
            current_node =  current_node->parent1();
        }
    }
}


void GraphBuilder::which_tip_is_below(){
    for ( auto it = nodes_.iterator(); it.good(); ++it){
        Node* current_node = (*it);
        if ( !current_node->isTip() ) continue;

        size_t tip_i = 0;
        for ( tip_i = 0 ; tip_i < current_node->tips_below.size(); tip_i++){
            if ( current_node->tips_below[tip_i] == 1 ) break;
        }

        while ( current_node->parent1() ){
            current_node->parent1()->tips_below[tip_i] = 1;
            current_node =  current_node->parent1();
        }
    }
}


void GraphBuilder::which_sample_is_below(){
    for ( auto it = nodes_.iterator(); it.good(); ++it){
        Node* current_node = (*it);
        if ( !current_node->isTip() ) continue;
        for ( size_t sample_i = 0 ; sample_i < current_node->samples_below.size(); sample_i++ ){
            Node * tmpNode = current_node;
            if ( current_node->samples_below[sample_i] == 1 ) {
                while ( current_node->parent1() ){
                    current_node->parent1()->samples_below[sample_i] = 1;
                    current_node =  current_node->parent1();
                }
            }
            current_node = tmpNode;
        }

    }
}


/*! \brief rewrite node content of nodes */
void GraphBuilder::rewrite_subTreeStr(){
    for ( size_t rank_i = 1; rank_i <= this->nodes_.back()->rank(); rank_i++){
        for ( auto it = nodes_.iterator(); it.good(); ++it){
            if ( (*it)->rank() != rank_i ) continue;

            (*it)->subTreeStr = ( (*it)->rank() == 1 ) ? (*it)->nodeName :
                                                         this->rewrite_internal_subTreeStr( (*it) );
        }
    }
}


string GraphBuilder::rewrite_internal_subTreeStr( Node * node ){
    string new_subTreeStr="(";
    for (size_t childIdx = 0; childIdx < node->child.size(); childIdx++ ){
        string bl1 = to_string ( node->child[childIdx]->edge1.bl() );
        if ( node->child[childIdx]->subTreeStr == node->child[childIdx]->nodeName ) {
            new_subTreeStr += node->child[childIdx]->nodeName+":" +  bl1;
        }
        else {
            bool new_hybrid_node = false;
            string brchlen_str2;
            for ( auto it = nodes_.iterator(); it.good(); ++it){
                if ( (*it) == node ) break;

                for ( size_t child_ii = 0; child_ii < (*it)->child.size(); child_ii++ ){
                    if ( (*it)->child[child_ii]->subTreeStr == node->child[childIdx]->subTreeStr){
                        new_hybrid_node = true;
                        brchlen_str2 = to_string(node->child[childIdx]->edge2.bl() );
                        break;
                    }
                }
                if (new_hybrid_node){
                    break;
                }
            }
            new_subTreeStr += new_hybrid_node ? node->child[childIdx]->nodeName+":" + brchlen_str2 :
                                                node->child[childIdx]->subTreeStr + node->child[childIdx]->nodeName+":" +  bl1;
        }
        if ( childIdx < node->child.size() - 1 ) new_subTreeStr += ",";
    }
    new_subTreeStr += ")";
    return new_subTreeStr;
}


void GraphBuilder::check_isUltrametric(){
    vector <Node*> remaining_node;//( this->nodes_.size(), 0 );
    for ( auto it = nodes_.iterator(); it.good(); ++it){
        remaining_node.push_back( (*it) );
    }
    size_t rank_i = 1;
    size_t remaining_node_i=0;
    while ( remaining_node.size() > 0 ){
        //int node_i = remaining_node[remaining_node_i];
        Node * current_node = remaining_node[remaining_node_i];
        if ( current_node->rank() == rank_i ){
            if (rank_i == 1) current_node->path_time.push_back(0.0);
            else{
                for (size_t child_i = 0; child_i < current_node->child.size(); child_i++ ){
                    double current_child_time = (current_node->child[child_i]->parent1() == current_node) ?
                                                current_node->child[child_i]->edge1.bl():
                                                current_node->child[child_i]->edge2.bl();
                    for ( size_t i = 0;i < current_node->child[child_i]->path_time.size(); i++){
                        current_node->path_time.push_back( current_child_time + current_node->child[child_i]->path_time[i] );
                    }
                }
            }
            remaining_node.erase( remaining_node.begin() + remaining_node_i );
        }
        else{
            remaining_node_i++;
        }

        if ( remaining_node_i == remaining_node.size()-1 ){
            rank_i++;
            remaining_node_i=0;
        }
    }

    for ( auto it = nodes_.iterator(); it.good(); ++it){
        for ( size_t i = 0; i < (*it)->path_time.size(); i++){
            if ( pow( ( (*it)->path_time[i] - (*it)->path_time[0] ) , 2 ) > 0.000001 ){
                this->is_ultrametric = false;
                break;
            }
        }
        //(*it)->set_height( (*it)->path_time[0] );
    }
}


string GraphBuilder::print_newick( Node * node ){
    string tree_str;
    if ( node->isTip() ) tree_str = node->nodeName ;
    else {
        tree_str = "(";
        for ( size_t i = 0 ; i < node->child.size() ; i++ ){
            tree_str += print_newick ( node->child[i] ) + ":" + to_string (node->child[i]->edge1.bl() );
            if ( i < node->child.size()-1 ) tree_str += ",";
        }
        tree_str += ")";
    }
    return tree_str;
}


/*! \brief Identify if its the start of the taxon name in a newick string, should be replaced by using (isalpha() || isdigit())  */
bool start_of_tax_name( string in_str, size_t i ){
    //bool start_bool = false;
    //if ( (in_str[i]!='(' && in_str[i-1]=='(') || (in_str[i-1]==',' && in_str[i]!='(') || ( (in_str[i-1]==')') && ( in_str[i]!=')' || in_str[i]!=':' || in_str[i]!=',' || in_str[i]!=';' ) ) ) {
        //start_bool=true;
    //}
    //return     start_bool;
    if      (  in_str[i-1] == '('  &&   in_str[i] != '(' ) return true;
    else if (  in_str[i-1] == ','  &&   in_str[i] != '(' ) return true;
    else if ( (in_str[i-1] == ')') && ( in_str[i] != ')' || in_str[i]!=':' || in_str[i]!=',' || in_str[i]!=';' ) ) return true;
    else return false;
}


size_t end_of_label_or_bl( string &in_str, size_t i ){
    for ( size_t j = i; j < in_str.size(); j++){
        if      ( in_str[j+1] == ',' )    return j;
        else if ( in_str[j+1] == ')' )    return j;
        else if ( in_str[j+1] == ':' )    return j;
        else if ( in_str[j+1] == ';' )    return j;
        else continue;
    }
}


void readNextStringto( string &readto, int& argc_i, int argc_, char * const* argv_ ){
    argc_i++;
    if (argc_i >= argc_) throw NotEnoughArg( argv_[argc_i-1] );
    readto = std::string(argv_[argc_i]);
    if ( readto[0] == '-' ) throw NotEnoughArg( argv_[argc_i-1] );
}


void GraphBuilder::removeOneChildInternalNode ( ){
    dout << "removeOneChildInternalNode";
    vector <Node*> toBeRemoved;
    for ( auto it = nodes_.iterator(); it.good(); ++it){
        if ( (*it)->child.size() != (size_t)1 || (*it)->isHybrid() ) {
            continue;
        }

        Node * removingFrom = (*it)->parent1();
        dout <<endl<< " Removing From parent " << removingFrom->nodeName << "( " << removingFrom->child.size() << " ) child: ";
        int removingChildIndex = -1;
        for ( size_t childIdx = 0; removingFrom->child.size(); childIdx++ ){
            if ( removingFrom->child[childIdx] == (*it) ){
                dout << (*it)->nodeName << endl;
                removingChildIndex = childIdx;
                break;
            }
        }

        assert (removingChildIndex != -1);
        removingFrom->child.erase(removingFrom->child.begin() + (size_t)removingChildIndex);
        (*it)->child[0]->set_parent1(NULL);
        removingFrom->add_child((*it)->child[0]);
        (*it)->child[0]->edge1.setLength( (*it)->child[0]->edge1.bl() + (*it)->edge1.bl());
        toBeRemoved.push_back((*it));
    }

    for ( size_t i = 0; i < toBeRemoved.size(); i++ ){
        nodes_.remove(toBeRemoved[i]);
    }
}


void GraphBuilder::removeZeroChildHybridNode(){
    dout << "removeZeroChildHybridNode";

    vector <Node*> toBeRemoved;
    for ( auto it = nodes_.iterator(); it.good(); ++it){
        if ( (*it)->child.size() != (size_t)0 || !(*it)->isHybrid() ) {
            continue;
        }

        Node * removingFrom = (*it)->parent1();
        dout <<endl<< " Removing From parent " << removingFrom->nodeName << "( " << removingFrom->child.size() << " ) child: ";
        int removingChildIndex = -1;
        for ( size_t childIdx = 0; removingFrom->child.size(); childIdx++ ){
            if ( removingFrom->child[childIdx] == (*it) ){
                dout << (*it)->nodeName << endl;
                removingChildIndex = childIdx;
                break;
            }
        }
        assert (removingChildIndex != -1);
        removingFrom->child.erase(removingFrom->child.begin() + (size_t)removingChildIndex);

        removingFrom = (*it)->parent2();
        dout <<endl<< " Removing From parent " << removingFrom->nodeName << "( " << removingFrom->child.size() << " ) child: ";
        removingChildIndex = -1;
        for ( size_t childIdx = 0; removingFrom->child.size(); childIdx++ ){
            if ( removingFrom->child[childIdx] == (*it) ){
                dout << (*it)->nodeName << endl;
                removingChildIndex = childIdx;
                break;
            }
        }
        assert (removingChildIndex != -1);
        removingFrom->child.erase(removingFrom->child.begin() + (size_t)removingChildIndex);

        toBeRemoved.push_back((*it));
    }

    for ( size_t i = 0; i < toBeRemoved.size(); i++ ){
        nodes_.remove(toBeRemoved[i]);
    }
}


//void GraphBuilder::mergeClade( valarray <size_t> clade ) {
    //vector <Node*> toBeRemoved;
    //for ( auto it = nodes_.iterator(); it.good(); ++it){
        //valarray <bool> comp = ( clade == (*it)->samples_below );
        //if ( comp.min() == true ){
            //for ( size_t child_i = 0; child_i < (*it)->child.size(); child_i ++){
                //Node * tmpNode = (*it)->child[child_i];
                //while (){
                    
                //}
            //}
            //break;
        //}
    //}
//}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*! \brief Remove interior nodes label of a string */
string remove_interior_label(string in_str/*!< input newick form string */){
    string out_str;
    out_str=in_str;

    size_t found_bracket=out_str.find(')');
    while ( found_bracket<out_str.size() ){
        if ( isalpha(out_str[found_bracket+1]) || isdigit(out_str[found_bracket+1]) ){
            size_t char_j = end_of_label_or_bl( out_str, found_bracket+1 );
            out_str.erase(out_str.begin()+found_bracket+1, out_str.begin()+char_j+1);
        }
        found_bracket = out_str.find( ")",found_bracket+1 );
    }
    return out_str;
}


//string rm_one_child_root(string in_str){
    //string out_str=tree_topo(in_str);
    //size_t first_bck_parent_idx=Parenthesis_balance_index_forwards(out_str,0);
    //size_t second_bck_parent_idx=Parenthesis_balance_index_forwards(out_str,1);
    //if ( (first_bck_parent_idx-1)== second_bck_parent_idx){
        //return in_str.substr(1,Parenthesis_balance_index_forwards(in_str,1))+in_str.substr(Parenthesis_balance_index_forwards(in_str,0)+1,in_str.size()-Parenthesis_balance_index_forwards(in_str,0)-1);
    //}
    //else{
        //return in_str;
    //}
//}


/////////////////////////////////////////// consider for removal

void GraphBuilder::rewrite_descendant(){    //check for coaleased tips(& sign in the tips)
    bool rewrite_descndnt=false;
    for ( auto it = nodes_.iterator(); it.good(); ++it){
        if ( (*it)->isTip() ){
            for (size_t i_str=0;i_str< (*it)->clade.size();i_str++){
                if ((*it)->clade[i_str]=='&'){
                    rewrite_descndnt=true;
                    break;
                }
            }
        }
        if (rewrite_descndnt){
            break;
        }
    }

    if ( !rewrite_descndnt ) return;

    //tax_name.clear();
    //int tax_name_start=0;
    //int tax_name_length=0;
    //for (size_t new_i_str=0;new_i_str<nodes_.back()->clade.size();new_i_str++){
        //tax_name_length++;
        //if (nodes_.back()->clade[new_i_str]=='&'){
            //tax_name_length--;
            //tax_name.push_back(nodes_.back()->clade.substr(tax_name_start,tax_name_length));
            //tax_name_start=new_i_str+1;
            //tax_name_length=0;
        //}
        //if (new_i_str==nodes_.back()->clade.size()-1){
            //tax_name.push_back(nodes_.back()->clade.substr(tax_name_start,tax_name_length));
        //}
    //}
    //sort(tax_name.begin(), tax_name.end());
    //descndnt.clear();

    //for ( auto it = nodes_.iterator(); it.good(); ++it){
        //vector <string> contained_tips;
        //valarray <int> re_initial_descndnt(0,tax_name.size());
        //int tax_name_start=0;
        //int tax_name_length=0;
        //for (size_t new_i_str=0;new_i_str<(*it)->clade.size();new_i_str++){
            //tax_name_length++;
            //if (nodes_.back()->clade[new_i_str]=='&'){
                //tax_name_length--;
                //contained_tips.push_back((*it)->clade.substr(tax_name_start,tax_name_length));
                //tax_name_start=new_i_str+1;
                //tax_name_length=0;
            //}
            //if (new_i_str== (*it)->clade.size()-1){
                //contained_tips.push_back( (*it)->clade.substr(tax_name_start,tax_name_length));
            //}
        //}
        //for (size_t tax_i=0;tax_i<tax_name.size();tax_i++){
            //for (size_t contained_tax_i=0;contained_tax_i<contained_tips.size();contained_tax_i++){
                //if (tax_name[tax_i]==contained_tips[contained_tax_i]){
                    ////descndnt[i][tax_i]=1;
                    //re_initial_descndnt[tax_i]=1;
                //}
            //}
        //}
        //descndnt.push_back(re_initial_descndnt);
    //}
    ////this->rewrite_node_clade();
    //this->init_node_clade();
}

#ifndef NDEBUG
bool GraphBuilder::print_all_node_dout(){

    if ( this->isNet() ) dout<<"           label  hybrid hyb_des non-tp parent1  height brchln1 parent2 brchln2 #child #dsndnt #id rank   edge   Clade "<<endl;
    else dout<<"       label non-tp   parent  height brchln #child #dsndnt #id rank edge   Clade "<<endl;
    for ( auto it = nodes_.iterator(); it.good(); ++it){
        //for (size_t j=0;j<descndnt[i].size();j++) dout<<setw(3)<<descndnt[i][j];

        assert( (*it)->print_dout( this->isNet() ) );

        dout<<"  ";

        //for (size_t j=0;j<samples_below[i].size();j++) dout<<samples_below[i][j];
        dout<<endl;
    }
    return true;
}
#endif
