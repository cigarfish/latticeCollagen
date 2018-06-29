////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File Name:  CSParameterContextTemporary.cpp                               //
//                                                                            //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                        //
//    Created:  2012-11-17 20:32:16                                           //
//                                                                            //
//  This file is part of the CellSys7 code.                                   //
//                                                                            //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig,                         //
//                   and INRIA, Paris-Rocquencourt.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "CSParameterContextTemporary.h"
#include "CSParameterTemporary.h"

#include <stack>


CSParameterContextTemporary::CSParameterContextTemporary( CSParameterContext * otherContext )
    : CSParameterContext( otherContext->name() )
{
    // cp children
    const std::vector<CSParameterContext *> & otherChildren =
        otherContext->getChildren();
    std::vector<CSParameterContext *>::const_iterator childrenIter =
        otherChildren.begin();

    while ( childrenIter != otherChildren.end() )
    {
        mChildren.push_back( new CSParameterContextTemporary( *childrenIter ) );
        ++childrenIter;
    }

    // cp Parameters
    const std::vector<CSParameter *> & otherParms = otherContext->getParameters();
    std::vector<CSParameter *>::const_iterator parmIt = otherParms.begin();

    while ( parmIt != otherParms.end() )
    {
        mParameters.push_back(  new CSParameterTemporary(*parmIt) );
        ++parmIt;
    }
}


CSParameterContextTemporary::CSParameterContextTemporary( const std::string & name )
    : CSParameterContext( name )
{}


CSParameterContextTemporary::~CSParameterContextTemporary()
{
    std::vector<CSParameter *> parms;
    std::vector<CSParameter *>::iterator parmIt;

    std::vector<CSParameterContext *> children;
    std::vector<CSParameterContext *>::iterator childIt;

    CSParameterContext * currentContext;

    // Put all sub-contexts on a stack to delete afterwards
    //  recursive deletion showed iterator invalidation (though I do not see why).

    // the stack on which all sub-contexts will be put:
    std::stack<CSParameterContext *> cStack;

    // an intermediate stack from which to pop a context when we process its subs:
    std::stack<CSParameterContext *> intermediateStack;

    intermediateStack.push( dynamic_cast<CSParameterContext *>(this) );

    childIt = mChildren.begin();

    while ( childIt != mChildren.end() )
    {
        intermediateStack.push( (*childIt) );
        ++childIt;
    }

    while ( ! intermediateStack.empty() )
    {
        currentContext = intermediateStack.top();
        intermediateStack.pop();

        cStack.push( currentContext );

        children = currentContext->getChildren();
        childIt = children.begin();

        while ( childIt != children.end() )
        {
            intermediateStack.push( *childIt );
            ++childIt;
        }

    }

    while ( !cStack.empty() )
    {
        currentContext = cStack.top();

        parms = currentContext->getParameters();
        parmIt = parms.begin();

        while ( parmIt != parms.end() )
        {
            delete (CSParameterTemporary *) (*parmIt);
            ++parmIt;
        }

        parms.clear();

        currentContext->getChildren().clear();

        if ( currentContext == this )
            return;

        delete currentContext;

        cStack.pop();
    }
}


CSParameterContextTemporary *
CSParameterContextTemporary::createFromXML( QXmlStreamReader * xml )
{
    // Todo:  Reimplement this for better memory management.
    return new CSParameterContextTemporary( CSParameterContext::createFromXML( xml ) );
}
