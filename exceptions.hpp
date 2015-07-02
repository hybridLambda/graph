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

#include <stdexcept>

#ifndef EXCEPTION
#define EXCEPTION

struct InvalidInput : public exception {
  string src;
  string reason;
  InvalidInput( ){
    this->src      = "";
    this->reason   = "";
  }
  InvalidInput( string str ){
    this->src      = str;
    this->reason   = "";
  }
  ~InvalidInput(){}
  virtual const char * what ( ) const throw () {
    string throwMsg = this->reason + this->src;
    return throwMsg.c_str();
  }
};


struct NotEnoughArg : public InvalidInput{
  NotEnoughArg( string str ):InvalidInput( str ){
    this->reason = "Not enough parameters when parsing option: ";
  }
  ~NotEnoughArg(){}
};


struct UnknowArg : public InvalidInput{
  UnknowArg( string str ):InvalidInput( str ){
    this->reason = "Unknow option: ";
  }
  ~UnknowArg(){}
};


struct GraphException : public InvalidInput{
  GraphException( string str ):InvalidInput( str ){}
  ~GraphException(){}
};


struct BranchLengthUnGiven : public GraphException{
  BranchLengthUnGiven( string str ):GraphException( str ){
    this->reason = "Branch length was not provided at node : ";
  }
  ~BranchLengthUnGiven(){}
};

struct ParenthesisNotBalanced : public GraphException{
  ParenthesisNotBalanced( string str ):GraphException( str ){
    this->reason = "Parenthesis not balanced: ";
  }
  ~ParenthesisNotBalanced(){}
};
#endif
