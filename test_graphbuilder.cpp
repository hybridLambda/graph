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

    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestGraphBuild );
