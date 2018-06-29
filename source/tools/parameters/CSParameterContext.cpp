///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSParameterContext.cpp                                               //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-03-14 17:04:44                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

#include "CSParameterContext.h"
#include "CSParameter.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stack>

#include <QtCore>

#include <QTreeView>
#include "../../gui/QCSParameterModel.h"
#include "../../gui/QCSParameterDelegate.h"


CSParameterContext::CSParameterContext( std::string name, CSParameterContext * parent )
    : CSParameterTreeItem( name, CSParameterTreeItem::Context, parent ),
      mpQtModel( NULL ),
      mpQtDelegate( NULL ),
      mDebug( false )
{
    mChildren.clear();
    mParameters.clear();

    // std::cout << "at the beginning..." << std::endl;
    // std::cout << "mChildren.size:  " << mChildren.size() << std::endl;
    // std::cout << "mParameters.size:  " << mParameters.size() << std::endl;
}


CSParameterContext::~CSParameterContext()
{
  clear();

  // unregister from parent context
  if ( mpParent )
    {
      std::vector<CSParameterContext *> siblings = ((CSParameterContext *)mpParent)->getChildren();
      std::vector<CSParameterContext *>::iterator foundMe =
        std::find( siblings.begin(), siblings.end(), this );
      if ( foundMe != siblings.end() )
      {
        if ( *foundMe == this  )
          ((CSParameterContext *)mpParent)->removeContext( this );
      }
      // else
      //   std::cout << "~CSParameterContext():  Context was not registered in parent!\n";
    }

  if ( mpQtModel )
      delete mpQtModel;
  if ( mpQtDelegate )
      delete mpQtDelegate;
}


void
CSParameterContext::rename( const std::string & newName )
{
  // does a context of this name already exist??
  // this should be looked after before calling this routine!
  mName = newName;
}


void
CSParameterContext::addContext(CSParameterContext * newChild)
{
    mChildren.push_back(newChild);
    // std::cout << "mChildren.size:  " << mChildren.size() << std::endl;
    newChild->setParent(this);
}


void
CSParameterContext::removeContext( CSParameterContext * obsoleteContext )
{
    std::vector<CSParameterContext *>::iterator found
      = std::find( mChildren.begin(), mChildren.end(), obsoleteContext );
    if( *found )
      mChildren.erase( found );
}

void
CSParameterContext::addParameter(CSParameter * newParm)
{
    if ( newParm )
    {
        mParameters.push_back(newParm);
        //std::cout << "mParameters.size:  " << mParameters.size() << std::endl;
        newParm->setParent(this);
    }
}


void
CSParameterContext::clearParameters()
{ mParameters.clear(); }


void
CSParameterContext::clear()
{
  // ouch, recursion..!
  mParameters.clear();
  std::vector<CSParameterContext *>::iterator subContext;
  for ( subContext = mChildren.begin();
        subContext != mChildren.end();
        ++subContext )
    {
      (*subContext)->clear();
    }

  mChildren.clear();
};

//!
// levelsDown == 0 means go down all the way
CSParameterContext *
CSParameterContext::findContext( const std::string & key, int levelsDown )
{
  // protruding breadth-first, putting non-fits on a stack (std::vector)
  std::vector<CSParameterContext *> stack;

  stack.push_back(this);

  int level = levelsDown;
  while ( stack.size() )
    {
      CSParameterContext * currentContext = stack.back();
      stack.pop_back();

      if ( currentContext->name() == key )
        return currentContext;

      if ( !levelsDown || level )
        if ( currentContext->hasChildren() )
          {
            std::vector<CSParameterContext *> children = currentContext->getChildren();
            stack.insert( stack.end(),
                          children.begin(),
                          children.end()
                          );
          }
      --level;
    }

  if ( mDebug )
      std::cout << "CSParameterContext::findContext:  Context \"" << key << "\" not found in Context \"" << name() << "\"\n";

  return NULL;
}


CSParameter *
CSParameterContext::findParameter( const std::string & key, const std::string & contextKey, int level )
{
  std::vector<CSParameterContext *> stack;
  CSParameterContext * inContext;

  if ( contextKey.size() == 0 )
    inContext = this;
  else
    if ( !(inContext = findContext( contextKey )) )
    {
      if ( mDebug )
        std::cout << "CSParameterContext::findParameter:  Context \""
                  << contextKey
                  << "\" not found in Context \""
                  << name() << "\"\n";
      return NULL;
    }

  int curLevel = 1;
  stack.push_back( inContext );

  while ( stack.size() )
    {
      CSParameterContext * curContext = stack.back();
      stack.pop_back();

      std::vector<CSParameter *> parms = curContext->getParameters();
      std::vector<CSParameter *>::const_iterator pIt;
      for ( pIt = parms.begin(); pIt != parms.end(); ++pIt )
        if ( (*pIt)->name() == key )
          return *pIt;

      if ( curLevel == level )
        break;

      if ( curContext->hasChildren() )
        {
          std::vector<CSParameterContext *> children = curContext->getChildren();
          stack.insert( stack.end(),
                        children.begin(),
                        children.end()
                        );
        }

      ++curLevel;
    }

  if ( mDebug )
    std::cout << "CSParameterContext::findParameter:  Parameter \""
              << key
              << "\" not found in Context \""
              << name() << "\"\n";

  return NULL;
}


std::vector<CSParameter *>
CSParameterContext::findParametersBySubstring( const std::string & key )
{
    std::vector<CSParameter *> returnList;
    std::string substring;

    substring = key;

    std::vector<CSParameterContext *> stack;

    int curLevel = 1;
    stack.push_back( this );

    while ( stack.size() )
    {
        CSParameterContext * curContext = stack.back();
        stack.pop_back();

        std::vector<CSParameter *> parms = curContext->getParameters();
        std::vector<CSParameter *>::const_iterator pIt;
        for ( pIt = parms.begin(); pIt != parms.end(); ++pIt ) {
            std::string paramName = (*pIt)->name();
            std::size_t found = paramName.find(substring);
            if ( found!=std::string::npos )
                returnList.push_back(*pIt);
        }

        if ( curContext->hasChildren() )
        {
            std::vector<CSParameterContext *> children = curContext->getChildren();
            stack.insert( stack.end(),
                          children.begin(),
                          children.end()
                );
        }

        ++curLevel;
    }

    return returnList;
}


std::vector<CSParameter *>
CSParameterContext::findParametersByDatatype( const CSParameter::DataType datatype )
{
    std::vector<CSParameter *> returnList;

    std::vector<CSParameterContext *> stack;

    int curLevel = 1;
    stack.push_back( this );

    while ( stack.size() )
    {
        CSParameterContext * curContext = stack.back();
        stack.pop_back();

        std::vector<CSParameter *> parms = curContext->getParameters();
        std::vector<CSParameter *>::const_iterator pIt;
        for ( pIt = parms.begin(); pIt != parms.end(); ++pIt ) {
            if ( (*pIt)->dataType() == datatype )
                returnList.push_back(*pIt);
        }

        if ( curContext->hasChildren() )
        {
            std::vector<CSParameterContext *> children = curContext->getChildren();
            stack.insert( stack.end(),
                          children.begin(),
                          children.end()
                );
        }

        ++curLevel;
    }

    return returnList;
}


#undef QT_NO_CAST_FROM_ASCII

void
CSParameterContext::writeXML( QXmlStreamWriter * xml, bool toplevel ) const
{
  if ( !toplevel )
  {
      xml->writeStartElement( "ParameterContext" );
      xml->writeAttribute( "name", name().c_str() );
  }

  std::vector<CSParameter *>::const_iterator parmsIt = mParameters.begin();
  while ( parmsIt != mParameters.end() )
    {
        (*parmsIt)->writeXML( xml );
        ++parmsIt;
    }
  // write all parameters of the context

  std::vector<CSParameterContext *>::const_iterator ctxtIt = mChildren.begin();
  while ( ctxtIt != mChildren.end() )
    {
      (*ctxtIt)->writeXML( xml );
      ++ctxtIt;
    }

  if ( !toplevel )
    xml->writeEndElement();
}


CSParameterContext *
CSParameterContext::createFromXML( QXmlStreamReader * xml )
{
    Q_ASSERT( xml->isStartElement() );

    std::string contextName = 
        xml->attributes().value("name").toString().toStdString();
    CSParameterContext * context = new CSParameterContext( contextName, NULL);

    while ( xml->readNextStartElement() )
    {
        if ( xml->name() == "Parameter" )
        {
            context->addParameter( CSParameter::createFromXML(xml) );
            xml->skipCurrentElement();
        }
        else if ( xml->name() == "ParameterContext" )
        {
            context->addContext( CSParameterContext::createFromXML(xml) );
            // std::cout << context->getChildren().size() << " subcontexts\n";
            // std::cout << context->getChildren()[0]->name() << std::endl;
        }
        else
        {
            std::cout << "Unknown Element:  " << xml->name().toString().toStdString() << std::endl;
            xml->skipCurrentElement();
        }
    }

    return context;
}


void
CSParameterContext::setupGUI( QTreeView * treeView )
{
    if ( ! treeView )
        return;

    if ( ! mpQtModel )
        mpQtModel = new QCSParameterModel( this );

    treeView->setModel( mpQtModel );

    if ( ! mpQtDelegate )
        mpQtDelegate = new QCSParameterDelegate();

    treeView->setItemDelegate( mpQtDelegate );

    QObject::connect( mpQtModel, SIGNAL( dataChanged(int) ),
                      treeView,  SLOT( resizeColumnToContents(int) ) );

    treeView->expandAll();

    // todo:  look how many columns there are:
    for ( int i = 0; i < 3; ++i )
        treeView->resizeColumnToContents( i );
}


void
CSParameterContext::dump( std::stringstream & output )
{
    std::stack<CSParameterContext *> conStack;

    conStack.push( this );

    while ( !conStack.empty() )
    {
        CSParameterContext * currentContext = conStack.top();
        conStack.pop();

        output << currentContext->name() << std::endl;


        const std::vector<CSParameter *> & parms = currentContext->getParameters();

        std::vector<CSParameter *>::const_iterator parmIt;

        for ( parmIt = parms.begin(); parmIt != parms.end(); ++parmIt )
        {
            output << (*parmIt)->name()
                   << "\t"
                   << (*parmIt)->dataString()
                   << std::endl;
        }


        const std::vector<CSParameterContext *> & contexts = currentContext->getChildren();

        std::vector<CSParameterContext *>::const_iterator conIt;

        for ( conIt = contexts.begin(); conIt != contexts.end(); ++conIt )
        {
            conStack.push( *conIt );
        }
    }
}


void
CSParameterContext::setVisible( bool visibility )
{
    std::stack<CSParameterContext *> conStack;

    conStack.push( this );

    while ( !conStack.empty() )
    {
        CSParameterContext * currentContext = conStack.top();
        conStack.pop();


        dynamic_cast<CSParameterTreeItem *>(currentContext)->setVisible( visibility );


        const std::vector<CSParameter *> & parms = currentContext->getParameters();
        std::vector<CSParameter *>::const_iterator parmIt;

        for ( parmIt = parms.begin(); parmIt != parms.end(); ++parmIt )
            (*parmIt)->setVisible( visibility );


        const std::vector<CSParameterContext *> & contexts = currentContext->getChildren();
        std::vector<CSParameterContext *>::const_iterator conIt;

        for ( conIt = contexts.begin(); conIt != contexts.end(); ++conIt )
        {
            conStack.push( *conIt );
        }
    }
}
