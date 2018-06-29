////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSRegistrar.h                                                 //
//                                                                            //
//     Author:  Andreas Buttenschoen <andreas@buttenschoen.ca>                //
//    Created:  2016-10-20                                                    //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <string>
#include "CSFactory.h"

template <class Interface, class T, typename... Args>
class CSRegistrar
{
public:
    CSRegistrar(const string& type)
    {
        //std::cout << "Registering " << type;
        auto instance = CSFactory<Interface, Args...>::Instance();
        //std::cout << " @ factory " << instance << std::endl;
        instance->template RegisterType<T>(type);
    }
};

// Factory method
template <class Interface, typename... Args>
auto create(const string& name, Args&&... args)
    -> typename CSFactory<Interface, Args...>::InterfacePointer
{
    auto instance = CSFactory<Interface, Args...>::Instance();
    return instance->Create(name, std::forward<Args>(args)...);
}


