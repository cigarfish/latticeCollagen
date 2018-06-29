///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  Lattice.h                                                            //
//                                                                                   //
//     Author:  Johannes Neitsch <johannes.neitsch@uni-leipzig.de>                   //
//    Created:  2012-10-11 16:20:03                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#ifndef LATTICE_H
#define LATTICE_H

#include <vector>

#include "../Elements/ModelElement.h"

//class ModelElement;

struct latticeContent{

  std::vector<ModelElement*> mpppvLattice;
  int index;

};

class Lattice
{

  public:
    Lattice(int xmin, int ymin, int zmin, int xl, int yl, int zl);

    void Add(ModelElement* modelElement);
    void Remove(ModelElement *modelElement);

    std::vector< ModelElement * > InteractWith(ModelElement *modelElement);
    bool check( ModelElement* modelElement);

    int mXminCut;
    int mYminCut;
    int mZminCut;

    int mXLenght;
    int mYLenght;
    int mZLenght;

    latticeContent ****mpppvLattice;

    int mXShift;
    int mYShift;
    int mZShift;

  private:


};


#endif //LATTICE_H
