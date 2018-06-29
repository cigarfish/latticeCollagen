
#include <cstring>

#include "Parameter.h"
#include "CSModelTools.h"

vector<string> tokenize(const string& str, const string& delimeters)
{
    vector<string> tokens;
    string::size_type fPos=0, lPos=0;
    do {
        fPos = str.find_first_not_of(delimeters, lPos);
        if (fPos == string::npos) break;
        lPos = str.find_first_of(delimeters, fPos);
        if (lPos == string::npos) lPos = str.length();
        tokens.push_back(str.substr(fPos, lPos - fPos));
    } while (lPos != str.length());

    return tokens;
}

XMLNode getXMLNode(const XMLNode XML_base, const string& path,
                   const string& attribute, const string& attributeNode)
{
    vector<string> tags = tokenize(path, "/");
    XMLNode xNode = XML_base;

    for (const auto tag : tags)
    {
        //std::cout << "tag=" << tag << std::endl;
        if (tag.empty()) continue;

        if (xNode.nChildNode(tag.c_str()) > 0)
            xNode = xNode.getChildNodeWithAttribute(tag.c_str(),
                                                   attribute.c_str(),
                                                   attributeNode.c_str());
      else
            return XMLNode::emptyXMLNode;

    }

    return xNode;
}

XMLNode getParameterNode(const XMLNode& XML_base, const string& name)
{
    return getXMLNode(XML_base, "Parameter", "name", name);
}

XMLNode getXMLNode(const XMLNode XML_base, const string& path)
{
    vector<string> tags=tokenize(path,"/");

    XMLNode xNode = XML_base;
    for (auto tag : tags)
    {
        if (tag.empty()) continue;
        int ntag = 0;
        string::size_type r_bracket, l_bracket;
        if ((r_bracket=tag.find_last_of("]")) != string::npos &&
            (l_bracket=tag.find_last_of("[")) != string::npos &&
            r_bracket > l_bracket)
        {
            stringstream tmp(tag.substr(l_bracket+1, r_bracket-1));
            tag.resize(l_bracket);
            tmp >> ntag;
        }
        if (xNode.nChildNode(tag.c_str()) > ntag)
            xNode = xNode.getChildNode(tag.c_str(), ntag);
        else
            return XMLNode::emptyXMLNode;
    }
    return xNode;
}

string strip_last_token(string& s, const string& del)
{
    string token;
    string::size_type pos = s.find_last_of(del);
    if (pos == string::npos)
    {
        token = s;
        s = "";
    }
    else
    {
        pos++;
        token = s.substr(pos);
        s.resize(pos-1);
    }
    return token;
}

//string remove_last_token(string& s, const string& del)
//{
//    string token;
//    string::size_type pos = s.find_last_of(del);
//    if (pos == string::npos)
//    {
//        token = s;
//        s = "";
//    }
//
//}

template <>
bool getXMLAttribute<string>(const XMLNode XML_base, string path, string& value,
                             bool verbose)
{
    string attribute;
    attribute = strip_last_token(path,"/");
    if (attribute.empty()) return false;

    if (verbose)
    {
        cout << "getXMLAttribute: seeking for: ";
        cout << "path:" << path << std::endl;
        if (path.empty()) cout << XML_base.getName();
        else cout << path;
        cout << "->" << attribute;
    }

    XMLNode xNode;
    if (path.length()>0)
    {
        xNode = getXMLNode(XML_base, path);
        if (verbose)
            std::cout << "1 Trying to access node with " << path << std::endl;
        //std::cout << "Node:" << xNode.getName() << std::endl;
    }
    else
    {
        xNode = XML_base;
    }

    XMLCSTR str_val;
    if (lower_case(attribute) == "text")
    {
        if (xNode.nText())
            str_val = xNode.getText();
        else
        {
            if (verbose)
                cout << " .. not found" << endl;
            return false;
        }
    }
    else if (xNode.isAttributeSet(attribute.c_str()))
        str_val = xNode.getAttribute(attribute.c_str());
    else
    {
        if (verbose)
            cout << " .. not found" << endl;
        return false;
    }
    value = string(str_val);
    if (path.length()==0)
        path = XML_base.getName();
    if (verbose)
        cout << ": " << value << endl;
    return true;
}

template<>
bool getXMLAttribute<bool>(const XMLNode XML_base, string path, bool& value,
                           bool verbose)
{
    string attribute;
    attribute = strip_last_token(path,"/");
    if (attribute.empty()) return false;

    if (verbose)
    {
        cout << "getXMLAttribute: seeking for: ";
        cout << "path:" << path << std::endl;
        if (path.length()==0) cout << XML_base.getName();
        else cout << path;
        cout << "->" << attribute << std::endl;
    }

    XMLNode xNode;
    if (!path.empty())
    {
        xNode = getXMLNode(XML_base, path);
        if (verbose)
            std::cout << "Trying to access node with " << path << std::endl;
        std::cout << "Node:" << xNode.getName() << std::endl;
    }
    else
    {
        xNode = XML_base;
    }

    XMLCSTR str_val;
    if (lower_case(attribute) == "text")
    {
        if (xNode.nText())
            str_val = xNode.getText();
        else
        {
            if (verbose)
                cout << " .. not found" << endl;
            return false;
        }
    }
    else if (xNode.isAttributeSet(attribute.c_str()))
        str_val = xNode.getAttribute(attribute.c_str());
    else
    {
        if (verbose)
            cout << " .. not found" << endl;
        return false;
    }

    // TODO how to avoid having to use strncmp? how to make sure null
    // terminations are uniform prior to comparision?
    str_val = lower_case(str_val).c_str();
    if ( ! strncmp(str_val, "1", 1) || ! strncmp(str_val, "true", 4) ||
         ! strncmp(str_val, "yes", 3))
    {
        value = true;
        if (verbose)
            cout << ": true" << endl;
        return true;
    }
    if ( ! strncmp(str_val, "0", 1) || ! strncmp(str_val, "false", 5) ||
         ! strncmp(str_val, "no", 2))
    {
        value = false;
        if (verbose)
            cout << ": false" << endl;
        return true;
    }

    cout << "Error converting values" << endl;
    return true;
}

