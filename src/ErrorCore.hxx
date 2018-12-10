//////////////////////////////////////////////////////////////
// $Id: EoaCore.hxx,v 1.5 2010/02/01 17:39:16 mcgrew Exp $
//
#ifndef EoaCore_HXX_SEEN
#define EoaCore_HXX_SEEN

#include <exception>
#include <string>

namespace COMET { class EoaCore; }

/// The root exception for all exceptions explicitly thrown by the oaCore
/// library.  All exceptions declared by this library should be derived from
/// this class, and not std::exception.  New exceptions can be derived from
/// EoaCore by creating a new class
///
/// \code
/// class EoaChild :public EoaCore {
/// public:
///     EoaChild() {AppendWhat("EoaChild");};
/// };
/// \endcode
///
/// But the prefered way to create derived exception classes is to use the
/// MYEXCEPTION macro.  This is used as follows
///
/// \code
/// MYEXCEPTION(EChild,EoaCore);
/// MYEXCEPTION(EGrandChild,EChild);
/// ...
/// try {
///   throw EChild();
/// }
/// catch (EoaCore& ex) {
///   std::cout << ex.what << std::endl;
/// }
/// \endcode
///
/// This can add a backtrace at the point of the exception.  The number of
/// backtrace frames to be shown is controlled with the
/// COMET::EoaCore::gBacktraceSymbols static variable.
///
class COMET::EoaCore :public std::exception {
    /// What exception generated this object.
    char fWhat[2048];

public:
    /// The number of backtrace symbols to add to the "what" string.
    static unsigned int gBacktraceSymbols;

    EoaCore();
    virtual ~EoaCore() throw() {}

    /// Inherited from exception to return the name of the exception as a null
    /// terminated string.
    const char* what(void) const throw() {return fWhat;}

    /// Used in constructors of classes which inherit EoaCore to add text to
    /// the What string.  Here is a comprehensive example of how to construct
    /// an exception derived from EoaCore:
    /// \verbatim
    /// class EoaChild : public EoaCore {
    /// public:
    ///   EoaChild() {AppendWhat("EoaChild");}
    /// };
    /// \endverbatim
    void AppendWhat(const char* child);
};

/// A macro to build an exception class __name that is derived from __parent.
/// The __parent class must be derived from EoaCore which provides the
/// AppendWhat method.  This macro builds an exception class that can be used:
/// \verbatim
/// MYEXCEPTION(EChild,EoaCore);
/// ...
/// try {
///   throw EChild();
/// }
/// catch (EoaCore& ex) {
///   std::cout << ex.what << std::endl;
/// }
/// \endverbatim
#ifndef MYEXCEPTION
#define MYEXCEPTION(__name,__parent)           \
    class __name : public __parent {            \
    public:                                     \
    __name() {AppendWhat(#__name);}             \
    }
#endif
#endif
