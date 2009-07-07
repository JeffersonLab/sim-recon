/*
 * Copyright 1999-2004 The Apache Software Foundation.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <xalanc/Include/PlatformDefinitions.hpp>

#include <cassert>

#if defined(XALAN_CLASSIC_IOSTREAMS)
#include <iostream.h>
#else
#include <iostream>
#endif

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>

#include <xalanc/PlatformSupport/XSLException.hpp>

//#include <xalanc/DOMSupport/XalanDocumentPrefixResolver.hpp>
#include <xalanc/PlatformSupport/PrefixResolver.hpp>
#include <xalanc/PlatformSupport/DOMStringHelper.hpp>

using namespace xercesc;
using namespace XALAN_CPP_NAMESPACE;

class myPrefixResolver : public PrefixResolver
{
 public:
    myPrefixResolver();
    myPrefixResolver(const myPrefixResolver &);
    ~myPrefixResolver();
    const XalanDOMString* getNamespaceForPrefix(
          const XalanDOMString &prefix) const;
    const XalanDOMString& getURI() const;

 private:
    XalanDOMString m_uri;
};

myPrefixResolver::myPrefixResolver()
: m_uri("http://www.gluex.org/hdds")
{
}

myPrefixResolver::myPrefixResolver(const myPrefixResolver &)
{
}

myPrefixResolver::~myPrefixResolver()
{
}

const XalanDOMString* myPrefixResolver::getNamespaceForPrefix(
                      const XalanDOMString &prefix) const
{
    return &m_uri;
}

const XalanDOMString& myPrefixResolver::getURI() const
{
    return m_uri;
}

#include <xalanc/XPath/XObject.hpp>
#include <xalanc/XPath/XPathEvaluator.hpp>

#include <xalanc/XalanSourceTree/XalanSourceTreeDOMSupport.hpp>
#include <xalanc/XalanSourceTree/XalanSourceTreeInit.hpp>
#include <xalanc/XalanSourceTree/XalanSourceTreeParserLiaison.hpp>

int main( int argc, char* argv[])
{
    int theResult = 0;

    if (argc != 4)
    {
        std::cerr << "Usage: SimpleXPathAPI XMLFilePath Context XPathExpression" << std::endl;
        theResult = -1;
    }
    else
    {
        try
        {
            XMLPlatformUtils::Initialize();
            XPathEvaluator::initialize();
            {
                // Initialize the XalanSourceTree subsystem...
                XalanSourceTreeInit        theSourceTreeInit;

                // We'll use these to parse the XML file.
                XalanSourceTreeDOMSupport        theDOMSupport;
                XalanSourceTreeParserLiaison    theLiaison(theDOMSupport);

                // Hook the two together...
                theDOMSupport.setParserLiaison(&theLiaison);

                const XalanDOMString    theFileName(argv[1]);

                // Create an input source that represents a local file...
                const LocalFileInputSource    theInputSource(theFileName.c_str());

                // Parse the document...
                XalanDocument* const    theDocument =
                        theLiaison.parseXMLStream(theInputSource);
                assert(theDocument != 0);

                myPrefixResolver    thePrefixResolver;
                XPathEvaluator    theEvaluator;

                // OK, let's find the context node...
                XalanNode* const    theContextNode =
                        theEvaluator.selectSingleNode(
                            theDOMSupport,
                            theDocument,
                            XalanDOMString(argv[2]).c_str(),
                            thePrefixResolver);

                if (theContextNode == 0)
                {
                    std::cerr << "Warning -- No nodes matched the location path \""
                         << argv[2]
                         << "\"."
                         << std::endl
                         << "Execution cannot continue..."
                         << std::endl
                         << std::endl;
                }
                else
                {
                    // OK, let's evaluate the expression...
                    const XObjectPtr    theResult(
                        theEvaluator.evaluate(
                                theDOMSupport,
                                theContextNode,
                                XalanDOMString(argv[3]).c_str(),
                                thePrefixResolver));

                    assert(theResult.null() == false);

                    std::cout << "The string value of the result is:"
                         << std::endl
                         << theResult->str()
                         << std::endl
                         << std::endl;
                }
            }

            XPathEvaluator::terminate();

            XMLPlatformUtils::Terminate();
        }
        catch(const XSLException&   theException)
        {
            std::cerr << "XSL exception: "
                 << theException.getMessage()
                 << std::endl;

            theResult = -1;
        }
        catch(...)
        {
            std::cerr << "Generic exception caught!" << std::endl;

            theResult = -1;
        }
    }

    return theResult;
}
