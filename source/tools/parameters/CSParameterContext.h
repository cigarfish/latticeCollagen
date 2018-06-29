///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//  File Name:  CSParameterContext.h                                                 //
//                                                                                   //
//     Author:  Tim Johann <tim.johann@uni-leipzig.de>                               //
//    Created:  2012-03-14 16:54:43                                                  //
//                                                                                   //
//  This file is part of the CellSys7 code.                                          //
//                                                                                   //
//  (C) The Drasdo Group, IZBI, Universitaet Leipzig, and INRIA, Paris-Rocquencourt. //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



#ifndef CS_PARAMETER_CONTEXT_H
#define CS_PARAMETER_CONTEXT_H

#include <string>
#include <vector>

#include "CSParameterTreeItem.h"
#include "CSParameter.h"
#include "CSParameterTemporary.h"

class QXmlStreamWriter;

class QTreeView;
class QCSParameterModel;
class QCSParameterDelegate;


class CSParameterContext : public CSParameterTreeItem
{
 public:

  //! Constructor.
  //! \param name   The name (and key) of the Context
  //! \param parent The parent context in which to register.
  CSParameterContext( const std::string name, CSParameterContext * parent =NULL );

  //! Destructor.
  //!  clears all sub-contexts and contained parameters and unregisters
  //!  from parent context.
  ~CSParameterContext();

  //! Return the parent context.
  //! Will return NULL if it is a root context.
  CSParameterContext * parent() const { return (CSParameterContext *) mpParent; };

  //! Set the parent context
  void setParent(CSParameterContext *parent)
  { CSParameterTreeItem::setParent( (CSParameterTreeItem *) parent ); };

  //! Add a child context explicitly.
  //! \param newChild The context to add.
  void addContext( CSParameterContext *newChild );

  //! Add a child context by name.
  //! \param contextName The name of the context to add.
  CSParameterContext * addContext( const std::string & contextName )
  {
      CSParameterContext *ctx = new CSParameterContext(contextName);
      addContext( ctx );
      return ctx;
  };

  //! Remove a context.
  //! \param obsoleteContext The context to remove from the list of children.
  //! This method will only remove the given obsoleteContext from the children,
  //!  leaving the obsoleteContext intact and in memory.
  void removeContext( CSParameterContext * obsoleteContext );

  //! Find a subContext by name
  //! \param key        The key (=name) of the context to find in the children or below.
  //! \param levelsDown How many levels to traverse.  A value of zero means all levels.
  CSParameterContext * findContext( const std::string & key, int levelsDown=0 );

  //! Return the list of sub-contexts.
  const std::vector<CSParameterContext *> & getChildren() const
    { return mChildren; };

  std::vector<CSParameterContext *> & getChildren()
    { return mChildren; };

  //! Add a parameter explicitly.
  //! \param newParm The parameter object to add.
  void addParameter( CSParameter * newParm );

  //! Add a parameter by specifying its data.
  //! \param parmName The name of the parameter (should be unique on its level within the context).
  //! \param type     The data type specified by an enum type from CSParameter::DataType.
  //! \param data     The pointer to the data container (cast to void *).
  //! \param unit     The string holding the data's unit.
  //! Creates a CSParameter with the given specifications.
  CSParameter * addParameter( const std::string parmName, CSParameter::DataType type,
                              void * data, const std::string unit )
  {
      CSParameter * ret = new CSParameter(parmName, type, data, unit, this);
      addParameter( ret );
      return ret;
  }

  //! Set the parameter's data, adding a new parameter if it does not exist in
  //! this context (excluding subcontexts).
  //! \param parmName The name of the parameter (should be unique on its level within the context).
  //! \param type     The data type specified by an enum type from CSParameter::DataType.
  //! \param data     The pointer to the data container (cast to void *).
  //! \param unit     The string holding the data's unit.
  //! Creates a CSParameter with the given specifications, if no parameter with
  //! the given parmName does not exists, otherwise copy the value of
  //!   void * data
  //! to the found parameter's mpData
  CSParameter * setParameter( const std::string parmName, CSParameter::DataType type,
                              void * data, const std::string unit )
  {
      if ( CSParameter * found = findParameter(parmName, 1) )
      {
          CSParameter parm(parmName, type, data, unit, this);
          *found = parm;
          return found;
      }
      else
      {
          CSParameter * parm = new CSParameter(parmName, type, data, unit, this);
          addParameter( parm );
          return parm;
      }
  }

  //! Clear the list of parameters.
  //! Clears the prameter list of this context.
  //! This does not clear the parameters of contexts contained in this context.
  //! It will also NOT delete() the dataPointers of the contexts.
  void clearParameters();

  //! Clear all data.
  //! Clears all parameters and contexts within this context after recursing through all sub-contexts. 
  void clear();
  

  //! Set the name of this context
  void rename( const std::string & );


  //! Return the list of parameters of this context.
  //! Does not traverse into sub-contexts.
  const std::vector<CSParameter *> & getParameters() const
    { return mParameters; };

  std::vector<CSParameter *> & getParameters()
    { return mParameters; };

  //! Find a parameter.
  //! \param key         Name of the parameter to find.
  //! \param contextKey  Name of the sub-context to search within.
  //! \param level       Maximum number of levels to descend and look for the key.
  //! Recursively searches for a parameter with the given key as name an returns
  //!  its pointer.  If level is not given or set to zero it will search through
  //!  all sub-contexts.
  //! Uses a breadth-first search through this context and its sub-contexts.
  CSParameter * findParameter( const std::string & key, const std::string & contextKey = "", int level=0 );

  //! Find a parameter.
  //! \param key         Name of the parameter to find.
  //! \param level       Maximum number of levels to descend and look for the key.
  //! Calls findParameter( key, "", level ).
  CSParameter * findParameter( const std::string & key, int level )
  { return findParameter( key, "", level ); };

  //! Find parameters containing a substring.
  //! \param key         Substring to be part of the parameter name.
  std::vector<CSParameter *> findParametersBySubstring( const std::string & key );

  //! Find parameters by datatype.
  //! \param key         Datatype of the parameters to find.
  std::vector<CSParameter *> findParametersByDatatype( const CSParameter::DataType datatype );

  //! Returns true if it has no parent context, false otherwise.
  bool isRoot() const { return (mpParent==NULL); };

  //! Return true, if this context has sub-contexts.
  bool hasChildren() const
  { return (mChildren.size()>0); };


  //! Write the XML description of this context.
  //! \param xmlStream The stream to write the XML description to.
  //! \param toplevel If the context is top level the start and end tags will not be written
  //! This method will write the description of the whole parameter tree
  //! recursing through all sub-contexts.
  void writeXML( QXmlStreamWriter *xmlStream, bool toplevel=false ) const;

  //! Create a CSParameterContext from an XML description.
  //! \param xmlStream The stream to read the description from.
  //! This method will create a CSParameterContext from the XML information in the xmlStream.
  //! It will create all CSParameters and sub-contexts, if the description contains them.
  static CSParameterContext * createFromXML( QXmlStreamReader * xmlStream );

  //! Setup the GUI for displaying the CSParameterContext's data.
  //! \param treeView The QTreeView for displaying the data.
  //! This routine takes care of creating the model and delegate for the given QTreeView,
  //! connecting signals and slots for updating the QTreeView's layout, when data has changed.
  void setupGUI( QTreeView * treeView );

  void dump( std::stringstream & output );

  void setVisible( bool visibility );

  void setDebug( bool debug =true ) { mDebug = debug; };


 protected:
  //! Sub-context list.
  std::vector<CSParameterContext *>  mChildren;
  //! Parameter list.
  std::vector<CSParameter *>       mParameters;

  QCSParameterModel * mpQtModel;
  QCSParameterDelegate * mpQtDelegate;

  bool mDebug;
};


#endif // CS_PARAMETER_CONTEXT_H
