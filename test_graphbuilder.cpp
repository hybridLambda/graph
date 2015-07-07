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


#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "graph.hpp"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestGraphBuild : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestGraphBuild );
  CPPUNIT_TEST ( testConstructor );
  CPPUNIT_TEST ( testRemoveOneChildInternalNode );
  CPPUNIT_TEST_SUITE_END();
 private:

 public:
  void setUp() {
        return;
    }

    void tearDown() {
        return;
    }

    void testConstructor(){
      //string treeStr = "((A:1,B:1):1,C:2);";

      string treeStr = "(C:2,(A:1,B:1):1);";
      CPPUNIT_ASSERT_NO_THROW ( GraphBuilder graph1(treeStr));
    }

    void testRemoveOneChildInternalNode(){
      string treeStr = "((C:1):1,(A:1,B:1):1);";
      //string treeStr = "(C:2,(A:1,B:1):1);";
      GraphBuilder graph1(treeStr);
      //CPPUNIT_ASSERT_NO_THROW( graph1.print() );
      CPPUNIT_ASSERT_NO_THROW( graph1.removeOneChildInternalNode() );
      //CPPUNIT_ASSERT_NO_THROW( graph1.print() );
      CPPUNIT_ASSERT_NO_THROW( graph1.rewrite_subTreeStr() );
      CPPUNIT_ASSERT_EQUAL ( string("((A:1.000000,B:1.000000)Int_2:1.000000,C:2.000000)"), graph1.subTreeStrAtRoot() );
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestGraphBuild );
