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


#include "node.hpp"

Node::Node ( size_t max_of_taxa,
             size_t max_of_tip,
             size_t max_of_sample, // number of tip
             string nodeName,
             string treeStr,
             double bl,
             bool tip ){
    this->init();
    this->nodeName = nodeName;
    this->subTreeStr = treeStr;
    this->edge1.setLength(bl);
    this->isTip_ = tip;
    // clade=" ";

    this->taxa_below = valarray < size_t > ( (size_t)0, max_of_taxa );
    this->tips_below = valarray < int > ( (int)0, max_of_tip );
    this->samples_below = valarray < size_t > ( (size_t)0, max_of_sample );
}

void Node::init(){
    this->parent1_ = NULL;
    this->parent2_ = NULL;
    this->previous_ = NULL;
    this->next_     = NULL;
    this->rank_     = 0;
    //this->height_ = 1.0/0.0;
    // prob_to_hybrid_left=1.0;
    this->visited_ = false;
    this->isBelowHybrid_ = false;
    this->num_descndnt=0;
    //this->num_descndnt_interior=0;
}


void Node::print( bool is_Net ){
    cout << setw(7) << this->nodeName;
    cout << setw(12) << (this);
    if ( is_Net ) cout << setw(6) << this->isHybrid();
    if ( is_Net ) cout << setw(8) << this->isBelowHybrid();
    cout << setw(5) << this->isTip();
    if ( this->parent1() != NULL ) cout << setw (11) << ( this->parent1() );
    else cout << "           ";
    //cout << setw (12) << this->height();
    cout << setw (12) << this->edge1.bl();
    if (is_Net){
        if ( this->parent2() != NULL) cout << setw (11) << ( this->parent2() );
        else cout << "           ";
        cout<<setw (12) << this->edge2.bl();
    }
    cout << setw (7) << this->child.size();
    cout << setw (8) << num_descndnt;
    cout << setw(4) << NumberOfInteriorNodesBelow();
    cout << setw(6) << this->rank() << "   ";

    //cout << setw(2)<<this->edge();
    //if ( is_Net ) cout << setw(3) << this->edge2();
    cout << "    " << this->clade;
    for ( size_t i = 0; i < this->tips_below.size(); i++ ){
        cout<<this->tips_below[i];
    }
    cout <<" ";
    for ( size_t i = 0; i < this->samples_below.size(); i++ ){
        cout<<this->samples_below[i];
    }
    cout <<" ";
    for ( size_t i = 0; i < this->taxa_below.size(); i++ ){
        cout<<this->taxa_below[i];
    }
    //cout<<endl;
}

/*! \brief Add child node to parent node */
void Node::add_child( Node *child_node /*! pointer to the child node*/){
    this->child.push_back(child_node);
    if ( child_node->parent1() != NULL ){
        child_node->set_parent2 ( this );
    }
    else child_node->set_parent1 ( this );
}


/*! \brief Rank network node from the bottom.
 *
 * Child node has lower rank than the parent node. Tip nodes have rank one, the root node has the highest rank
 */
void Node::CalculateRank(){
    if ( this->isTip() ) {
        this->rank_ = 1;
        return;
    }
    else {
        size_t child_max_rank = 0;
        for ( size_t ith_child = 0; ith_child < this->child.size(); ith_child++ ){
            this->child[ith_child]->CalculateRank();
            child_max_rank = max( child_max_rank, this->child[ith_child]->rank() );
        }
        this->rank_ = child_max_rank + 1;
        return;
    }
}




///*! \brief Label a node if its a descendant of a hybrid node */
void Node::findWhoIsBelowHybrid(){
    if ( this->isTip() )
        return;

    for ( size_t childIndex = 0; childIndex < this->child.size(); childIndex++){
        if ( this->isHybrid() || this->isBelowHybrid() )
            this->child[childIndex]->setIsBelowHybrid( true );

        this->child[childIndex]->findWhoIsBelowHybrid();
    }
}



#ifndef NDEBUG

bool Node::print_dout( bool is_Net ){
    dout << setw(12) << this;
    if ( is_Net ) dout << setw(6) << this->isHybrid();
    if ( is_Net ) dout << setw(8) << this->isBelowHybrid();
        dout << setw(5) << this->isTip();
    if ( this->parent1() != NULL ) dout << setw (11) << ( this->parent1() ); //if (this->parent1) dout << setw (11) << (this->parent1_());
    else dout << "           ";
    //dout << setw (6) << this->height();
    //dout << setw (12) << this->brchlen1();
    //if (is_Net){
        //if ( this->parent2() != NULL ) dout << setw (11) << ( this->parent2()->nodeName );
        //else dout << "           ";
        //dout<<setw (12) << this->brchlen2();
    //}
    dout << setw (7) << this->child.size();
    dout << setw (8) << num_descndnt;
    dout << setw(4) << NumberOfInteriorNodesBelow();
    dout << setw(6) << this->rank() << "   ";
    //for (size_t i=0;i<descndnt.size();i++){
        //dout<<setw (1)<<descndnt[i];
    //}
    //dout << setw(2)<<this->edge();
    //if ( is_Net ) dout << setw(3) << this->edge2();
    dout << "    " << this->clade;
    //dout<<endl;
    return true;
}

#endif
