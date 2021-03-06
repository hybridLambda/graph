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
#include <stdexcept>
#include <cfloat>
#include <cassert>
//#include <iostream>

#ifndef NODECONTAINER
#define NODECONTAINER

//class NodeContainer;

class NodeIterator;
class ConstNodeIterator;

class NodeContainer {
    friend class GraphBuilder;
    friend class NodeIterator;
    friend class ConstNodeIterator;
    friend class CoalGT;
    friend class CoalST;
    friend class CoalSN;
    friend class TmpSN;
    friend class Figure;
    //public:
    NodeContainer();
    ~NodeContainer() { this->clear(); };

    NodeContainer& operator=( NodeContainer nc ) {
        swap(*this, nc);
        return(*this);
    };

    NodeContainer( const NodeContainer &nc );

    NodeIterator iterator();
    NodeIterator iterator( Node* node );
    ConstNodeIterator iterator() const;
    ConstNodeIterator iterator(Node* node) const;

    //private:
    friend void swap( NodeContainer& a, NodeContainer& b );

    size_t size_;
    Node* first_node_;
    Node* last_node_;

    void add( Node* node );
    void remove( Node *node );
    //void move( Node *node, const double new_height );
    void clear();

    Node* at(size_t nr) const;
    Node const* get(size_t nr) const { return at(nr); };

    Node* first() const { return first_node_; };
    Node* back() const { return last_node_; };

    size_t size() const { return size_; };

    void set_first( Node* node) { this->first_node_ = node; }
    void set_last ( Node* node) { this->last_node_  = node; }
    void add_before(Node* add, Node* before); // Use this when adding new nodes.

};

class NodeIterator {
    Node* current_node_;
  public:

    NodeIterator() { this->current_node_ = NULL; };
    NodeIterator( NodeContainer& nc) { this->current_node_ = nc.first(); };
    NodeIterator( Node* node) { this->current_node_ = node; };
    ~NodeIterator() {};

    Node* operator*() {
      if ( this->current_node_ == NULL ) throw std::out_of_range( "Node iterator out of range" );
      return this->current_node_;
    }

    Node* operator++() {
      if ( this->current_node_ == NULL ) throw std::out_of_range( "Node iterator out of range" );
      this->current_node_ = ( this->current_node_->is_last() ) ? NULL : this->current_node_->next();
      return this->current_node_;
    }

    Node* operator--() {
      if ( this->current_node_ == NULL ) throw std::out_of_range( "Node iterator out of range" );
      this->current_node_ = ( this->current_node_->is_first() ) ? NULL : this->current_node_->previous();
      return this->current_node_;
    }

    bool good() const { return ( this->current_node_ != NULL ); }

    //double height() const { return ( this->good() ? this->current_node_->height() : DBL_MAX ); }
};

class ConstNodeIterator {
    Node const* current_node_;
    size_t node_index_;
  public:
    size_t node_index() const { return this->node_index_; }

    ConstNodeIterator() { current_node_ = NULL; };
    ConstNodeIterator( const NodeContainer& nc) { current_node_ = nc.first(); };
    ConstNodeIterator( Node const* node) { current_node_ = node; };
    ~ConstNodeIterator() {};

    Node const* operator*() {
      if ( current_node_ == NULL ) throw std::out_of_range( "Node iterator out of range" );
      return current_node_;
    }

    Node const* operator++() {
      if ( current_node_ == NULL ) throw std::out_of_range( "Node iterator out of range" );
      current_node_ = ( current_node_->is_last() ) ? NULL : current_node_->next();
      return current_node_;
    }

    Node const* operator--() {
      if ( current_node_ == NULL ) throw std::out_of_range( "Node iterator out of range" );
      current_node_ = ( current_node_->is_first() ) ? NULL : current_node_->previous();
      return current_node_;
    }

    bool good() const { return ( current_node_ != NULL ); }

    //double height() const { return ( this->good() ? this->current_node_->height() : DBL_MAX ); }
};



inline NodeIterator NodeContainer::iterator() { return NodeIterator(*this); }
inline NodeIterator NodeContainer::iterator(Node* node) { return NodeIterator(node); }
inline ConstNodeIterator NodeContainer::iterator() const { return ConstNodeIterator(*this); }
inline ConstNodeIterator NodeContainer::iterator(Node* node) const { return ConstNodeIterator(node); }
#endif
