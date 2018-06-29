////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSFactory.h                                                   //
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
#include <memory>
#include <string>
#include <unordered_map>
#include <functional>
#include <stdexcept>

#include "../utils/make_unique.h"
#include "type_info.h"
#include "Types.h"

using std::string;
using std::function;
using std::unique_ptr;
using std::shared_ptr;
using std::unordered_map;

template <class Interface, typename... Args>
class CSFactory
{
public:
    using InterfacePointer = std::unique_ptr<Interface>;

    static CSFactory * Instance()
    {
        static CSFactory factory;
        return &factory;
    }

    template <class T>
    void RegisterType(const string& type)
    {
        auto it(mCreatorMap.find(type));

        //if (it != mCreatorMap.end())
        //{
        //    std::cout << "Type '" << type << "' already has registered a factory"
        //        << std::endl;
        //    std::cerr << "FACTORY AT " << this << " for type "
        //        << type_name<Args...>() << std::endl;
        //}

        mCreatorMap[type].reset(new Creator<T>);
    }

    InterfacePointer Create(const string& type, Args&&... args)
    {
        auto it(mCreatorMap.find(type));

        //std::cerr << "FACTORY AT " << this << " for type "
        //    << type_name<Args...>() << std::endl;

        if (it != mCreatorMap.end())
        {
            return (it->second)->create(std::forward<Args>(args)...);
        }


        throw std::logic_error("Cannot find factory for type '" + type + "'");
    }

private:
    CSFactory() {}

    struct ICreator
    {
        virtual InterfacePointer create(Args&&... args) = 0;
    };

    template <class T>
    struct Creator : public ICreator
    {
        virtual InterfacePointer create(Args&&... args) override
        {
            return std::make_unique<T>(std::forward<Args>(args)...);
        }
    };

    // map to the factories
    unordered_map<string, std::unique_ptr<ICreator>> mCreatorMap;
};
