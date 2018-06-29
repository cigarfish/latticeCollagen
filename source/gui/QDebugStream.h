// QDebugStream.h
// derived from:

// qdebugstream.h          20.7.10
// version 1.1             4.11.11

// version 1.1: patch ("size_t pos") to fix 64-bit error

// from Jochen Ulrich 6.6.2005 posting to qt-interest at trolltech.com
// "Redirecting std::cout to QTextEdit / Using QTextEdit as a Log Window"
// http://lists.trolltech.com/qt-interest/2005-06/thread00166-0.html

// example usage from same posting:

//    #include "qdebugstream.h"
//    #include "qtextedit.h"
//
//    void main( )
//    {
//       [...]
//       QTexEdit* myTextEdit = new QTextEdit(this, "myTextEdit");
//       myTextEdit->setTextFormat(Qt::LogText);
//
//       QDebugStream qout(std::cout, myTextEdit);
//       std::cout << "Send this to the Text Edit!" << endl;
//
//       [...]
//    }

// In this derivative, we want to duplex the input to the ostream to the
// QTextEdit, rather than just redirect.  I.e. the output goes to both the
// ostreams streambuf and the QTextEdit.
// This is not safe, as in the case of an ofstream, a close() will not close
// the file.  This is due to the kidnapping of the ostream's streambuf

//################
//# qdebugstream.h  #
//################

#ifndef Q_DEBUG_STREAM_H
#define Q_DEBUG_STREAM_H

#include <iostream>
#include <streambuf>
#include <string>

#include <QTextEdit>

class QDebugStream : public std::basic_streambuf<char>
{
public:
 QDebugStream(std::ostream &stream, QTextEdit* text_edit) : m_stream(stream), orig_ostream(stream.rdbuf())
 {
  log_window = text_edit;
  m_old_buf = stream.rdbuf();
  stream.rdbuf(this);
 }
 ~QDebugStream()
 {
  // output anything that is left
   if (!m_string.empty()) {
     log_window->append(m_string.c_str());
     orig_ostream << m_string;
     orig_ostream.flush();
   }
   orig_ostream.rdbuf(this);
   m_stream.rdbuf(m_old_buf);
 }

protected:
 virtual int_type overflow(int_type v)
 {
  if (v == '\n')
  {
   log_window->append(m_string.c_str());
   orig_ostream << m_string << std::endl;
   m_string.erase(m_string.begin(), m_string.end());
  }
  else
   m_string += v;

  return v;
 }

 virtual std::streamsize xsputn(const char *p, std::streamsize n)
 {
  m_string.append(p, p + n);

  size_t pos = 0;  // patch to avoid 64-bit conversion error
  while (pos != std::string::npos)
  {
   pos = m_string.find('\n');
   if (pos != std::string::npos)
   {
    std::string tmp(m_string.begin(), m_string.begin() + pos);
    log_window->append(tmp.c_str());
    orig_ostream << tmp << std::endl;
    m_string.erase(m_string.begin(), m_string.begin() + pos + 1);
   }
  }

  return n;
 }

private:
 std::ostream &m_stream;
 std::ostream orig_ostream;
 std::streambuf *m_old_buf;
 std::string m_string;
 QTextEdit* log_window;
};

#endif
