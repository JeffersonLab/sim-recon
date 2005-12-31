<?xml version="1.0" encoding="UTF-8"?>
<xsl:transform xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
               xmlns:xs="http://www.w3.org/2001/XMLSchema"
               xmlns="http://www.gluex.org/hddm"
               exclude-result-prefixes="xs xsl"
               version="1.0">

              <!-- xmlns:hddm="http://www.gluex.org/hddm"-->
<xsl:output method="xml" version="1.0" encoding="iso-8859-1" indent="no"/>

<!-- Top-level processing: make sure document tag is xs:schema -->

<xsl:template match="/">
  <xsl:for-each select="*">
    <xsl:if test="local-name() != 'schema'">
      <xsl:message terminate="yes">
  Format error: expected schema, found <xsl:value-of select="name()"/>
      </xsl:message>
    </xsl:if>
    <xsl:if test="namespace-uri() != 'http://www.w3.org/2001/XMLSchema'">
      <xsl:message terminate="yes">
  Format error: wrong namespace for schema:
       expected "http://www.w3.org/2001/XMLSchema"
       found <xsl:value-of select="namespace-uri()"/>
      </xsl:message>
    </xsl:if>
  </xsl:for-each>
  <xsl:apply-templates select="/xs:schema"/>
</xsl:template>

<!-- check that the document tag is hddm:HDDM -->

<xsl:template match="/xs:schema">
  <xsl:if test="count(xs:element[@name='HDDM']) = 0">
    <xsl:message terminate="yes">
  Format error: schema does not describe HDDM document type
    </xsl:message>
  </xsl:if>
  <xsl:if test="count(xs:element[@name='HDDM']) > 1">
    <xsl:message terminate="yes">
  Format error: element HDDM should only appear once in schema
    </xsl:message>
  </xsl:if>
  <xsl:if test="@targetNamespace != 'http://www.gluex.org/hddm'">
    <xsl:message terminate="yes">
  Format error: target namespace should be "http://www.gluex.org/hddm"
    </xsl:message>
  </xsl:if>
  <xsl:for-each select="/xs:schema/xs:element[@name='HDDM']">
    <xsl:call-template name="HDDM"/>
  </xsl:for-each>
</xsl:template>

<!-- check HDDM tag attributes and unroll template -->

<xsl:template name="HDDM">
  <xsl:variable name="class">
    <xsl:value-of select="xs:complexType/xs:attribute[@name='class']/@fixed"/>
  </xsl:variable>
  <xsl:if test="$class = ''">
    <xsl:message terminate="yes">
  Format error: unable to find HDDM "class" attribute
    </xsl:message>
  </xsl:if>
  <xsl:variable name="version">
    <xsl:value-of select="xs:complexType/xs:attribute[@name='version']/@fixed"/>
  </xsl:variable>
  <xsl:if test="$version = ''">
    <xsl:message terminate="yes">
  Format error: unable to find HDDM "version" attribute
    </xsl:message>
  </xsl:if>
  <xsl:for-each select="xs:complexType/xs:attribute">
    <xsl:choose>
      <xsl:when test="@name = 'class'"/>
      <xsl:when test="@name = 'version'"/>
      <xsl:otherwise>
        <xsl:message terminate="yes">
  Format error: schema defines unknown attribute "<xsl:value-of select="@name"/>"
  for HDDM element
        </xsl:message>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:for-each>
  <xsl:call-template name="element"/>
</xsl:template>

<!-- generate user tags -->

<xsl:template name="element" match="xs:element[@name]">
  <xsl:param name="minOccurs" select="1"/>
  <xsl:param name="maxOccurs" select="1"/>
  <xsl:for-each select="*">
    <xsl:choose>
      <xsl:when test="name() = 'xs:complexType' and
                      namespace-uri() = 'http://www.w3.org/2001/XMLSchema'"/>
      <xsl:otherwise>
        <xsl:message terminate="yes">
  Format error: elements of type "<xsl:value-of select="name()"/>"
                from namespace "<xsl:value-of select="namespace-uri()"/>"
  are not currently supported
        </xsl:message>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:for-each>
  <xsl:for-each select="@*">
    <xsl:choose>
      <xsl:when test="name() = 'name'"/>
      <xsl:when test="name() = 'id'"/>
      <xsl:when test="name() = 'nillable'"/>
      <xsl:when test="name() = 'minOccurs'">
        <xsl:message terminate="yes">
  Format error: element references not allowed in a top-level declaration
        </xsl:message>
      </xsl:when>
      <xsl:when test="name() = 'maxOccurs'">
        <xsl:message terminate="yes">
  Format error: element references not allowed in a top-level declaration
        </xsl:message>
      </xsl:when>
      <xsl:when test="name() = 'ref'">
        <xsl:message terminate="yes">
  Format error: element references not allowed in a top-level declaration
        </xsl:message>
      </xsl:when>
      <xsl:when test="name() = 'substitutionGroup'">
        <xsl:message terminate="yes">
  Format error: substitution groups are not allowed by HDDM
        </xsl:message>
      </xsl:when>
      <xsl:when test="name() = 'fixed'">
        <xsl:message terminate="yes">
  Format error: fixed elements (as opposed to fixed attributes)
  are not allowed by HDDM
        </xsl:message>
      </xsl:when>
      <xsl:when test="name() = 'default'">
        <xsl:message terminate="yes">
  Format error: defaulted elements are not allowed by HDDM
        </xsl:message>
      </xsl:when>
      <xsl:otherwise>
        <xsl:message terminate="yes">
  Format error: element attribute "<xsl:value-of select="name()"/>"
  is not currently supported in the element declaration
        </xsl:message>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:for-each>
  <xsl:element name="{@name}">
    <xsl:if test="$minOccurs != 1">
      <xsl:attribute name="minOccurs">
        <xsl:value-of select="$minOccurs"/>
      </xsl:attribute>
    </xsl:if>
    <xsl:if test="$maxOccurs != 1">
      <xsl:attribute name="maxOccurs">
        <xsl:value-of select="$maxOccurs"/>
      </xsl:attribute>
    </xsl:if>
    <xsl:choose>
      <xsl:when test="count(./xs:complexType) > 0">
        <xsl:apply-templates/>
      </xsl:when>
      <xsl:when test="count(@type) > 0">
        <xsl:variable name="userType">
          <xsl:value-of select="@type"/>
        </xsl:variable>
        <xsl:if test="count(ancestor::node()/xs:complexType[@name=$userType]) = 1">
          <xsl:message terminate="yes">
  Format error: user-defined type "<xsl:value-of select="$userType"/>"
                assigned to element "<xsl:value-of select="name()"/>"
  must be a named complexType to fit the HDDM specification
          </xsl:message>
        </xsl:if>
        <xsl:apply-templates select="ancestor::node()/xs:complexType[@name=$userType]"/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:message terminate="yes">
  Format error: element "<xsl:value-of select="name()"/>"
  has a type declaration that does not fit the HDDM specification
        </xsl:message>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:element>
</xsl:template>

<xsl:template match="xs:complexType">
  <xsl:for-each select="*">
    <xsl:choose>
      <xsl:when test="name() = 'xs:sequence' and
                      namespace-uri() = 'http://www.w3.org/2001/XMLSchema'"/>
      <xsl:when test="name() = 'xs:attribute' and
                      namespace-uri() = 'http://www.w3.org/2001/XMLSchema'"/>
      <xsl:when test="name() = 'xs:attributeGroup' and
                      namespace-uri() = 'http://www.w3.org/2001/XMLSchema'"/>
      <xsl:otherwise>
        <xsl:message terminate="yes">
  Format error: contents of type "<xsl:value-of select="name()"/>"
                from namespace "<xsl:value-of select="namespace-uri()"/>"
  are not currently supported
        </xsl:message>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:for-each>
  <xsl:apply-templates select="xs:attribute"/>
  <xsl:apply-templates select="xs:attributeGroup"/>
  <xsl:apply-templates select="xs:sequence"/>
</xsl:template>

<xsl:template match="xs:sequence">
  <xsl:for-each select="*">
    <xsl:choose>
      <xsl:when test="name() = 'xs:element' and
                      namespace-uri() = 'http://www.w3.org/2001/XMLSchema'"/>
      <xsl:otherwise>
        <xsl:message terminate="yes">
  Format error: sequence contains "<xsl:value-of select="name()"/>"
                from namespace "<xsl:value-of select="namespace-uri()"/>"
  not currently supported
        </xsl:message>
      </xsl:otherwise>
    </xsl:choose>
    <xsl:if test="count(@ref) = 0">
      <xsl:message terminate="yes">
  Format error: sequence contains anonymous "<xsl:value-of select="name()"/>"
  all elements must be named and declared at the top level of the schema
      </xsl:message>
    </xsl:if>
  </xsl:for-each>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template name="elementRef" match="xs:element[@ref]">
  <xsl:for-each select="@*">
    <xsl:choose>
      <xsl:when test="name() = 'ref'"/>
      <xsl:when test="name() = 'minOccurs'"/>
      <xsl:when test="name() = 'maxOccurs'"/>
      <xsl:when test="name() = 'id'"/>
      <xsl:when test="name() = 'nillable'"/>
      <xsl:when test="name() = 'abstract'">
        <xsl:message terminate="yes">
  Format error: abstract element not allowed in an element reference
        </xsl:message>
      </xsl:when>
      <xsl:when test="name() = 'type'">
        <xsl:message terminate="yes">
  Format error: type declaration not allowed in an element reference
        </xsl:message>
      </xsl:when>
      <xsl:when test="name() = 'substitutionGroup'">
                      namespace-uri() = 'http://www.w3.org/2001/XMLSchema'">
        <xsl:message terminate="yes">
  Format error: substitution groups are not allowed by HDDM
        </xsl:message>
      </xsl:when>
      <xsl:when test="name() = 'fixed'">
        <xsl:message terminate="yes">
  Format error: fixed elements (as opposed to fixed attributes)
  are not allowed by HDDM
        </xsl:message>
      </xsl:when>
      <xsl:when test="name() = 'default'">
        <xsl:message terminate="yes">
  Format error: defaulted elements are not allowed by HDDM
        </xsl:message>
      </xsl:when>
      <xsl:otherwise>
        <xsl:message terminate="yes">
  Format error: element attribute "<xsl:value-of select="name()"/>"
  is not currently supported in an element reference
        </xsl:message>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:for-each>
  <xsl:variable name="userTag">
    <xsl:value-of select="substring-after(@ref,':')"/>
  </xsl:variable>
  <xsl:if test="count(/xs:schema/xs:element[@name=$userTag]) != 1">
    <xsl:message terminate="yes">
  Format error: element "<xsl:value-of select="$userTag"/>"
  must be defined once at the top level of the document
    </xsl:message>
  </xsl:if>
  <xsl:apply-templates select="/xs:schema/xs:element[@name=$userTag]">
    <xsl:with-param name="minOccurs" select="@minOccurs"/>
    <xsl:with-param name="maxOccurs" select="@maxOccurs"/>
  </xsl:apply-templates>
</xsl:template>

<xsl:template match="xs:attribute">
  <xsl:if test="count(@use) > 0 and @use != 'required'">
    <xsl:message terminate="no">
  Warning: use="<xsl:value-of select="@use"/>" ignored
  for attribute "<xsl:value-of select="name(current())"/>
    </xsl:message>
  </xsl:if>
  <xsl:choose>
    <xsl:when test="count(@ref) > 0">
      <xsl:for-each select="@*">
        <xsl:choose>
          <xsl:when test="name() = 'id'"/>
          <xsl:when test="name() = 'use'"/>
          <xsl:otherwise>
            <xsl:message terminate="yes">
  Format error: attribute reference "<xsl:value-of select="name(current())"/>"
  should not try to define its own "<xsl:value-of select="name()"/>"
            </xsl:message>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:for-each>
      <xsl:variable name="userAtt">
        <xsl:value-of select="@ref"/>
      </xsl:variable>
      <xsl:apply-templates select="ancestor::node()/xs:attribute[@name=$userAtt]"/>
    </xsl:when>
    <xsl:when test="count(@name) > 0">
      <xsl:for-each select="@*">
        <xsl:choose>
          <xsl:when test="name() = 'name'"/>
          <xsl:when test="name() = 'fixed'"/>
          <xsl:when test="name() = 'id'"/>
          <xsl:when test="name() = 'type'"/>
          <xsl:when test="name() = 'use'"/>
          <xsl:otherwise>
            <xsl:message terminate="yes">
  Format error: attribute descriptor "<xsl:value-of select="name()"/>"
  is not currently supported
            </xsl:message>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:for-each>
      <xsl:choose>
        <xsl:when test="count(@fixed)">
          <xsl:attribute name="{@name}">
            <xsl:value-of select="@fixed"/>
          </xsl:attribute>
        </xsl:when>
        <xsl:when test="count(@type)">
          <xsl:variable name="thisType">
            <xsl:value-of select="@type"/>
          </xsl:variable>
          <xsl:variable name="localType">
            <xsl:value-of select="substring-after($thisType,':')"/>
          </xsl:variable>
          <xsl:choose>
            <xsl:when test="$thisType = 'hddm:int'">
              <xsl:attribute name="{@name}">int</xsl:attribute>
            </xsl:when>
            <xsl:when test="$thisType = 'hddm:long'">
              <xsl:attribute name="{@name}">long</xsl:attribute>
            </xsl:when>
            <xsl:when test="$thisType = 'hddm:float'">
              <xsl:attribute name="{@name}">float</xsl:attribute>
            </xsl:when>
            <xsl:when test="$thisType = 'hddm:double'">
              <xsl:attribute name="{@name}">double</xsl:attribute>
            </xsl:when>
            <xsl:when test="$thisType = 'hddm:Particle_t'">
              <xsl:attribute name="{@name}">Particle_t</xsl:attribute>
            </xsl:when>
            <xsl:when test="$thisType = 'hddm:string'">
              <xsl:attribute name="{@name}">string</xsl:attribute>
            </xsl:when>
            <xsl:when test="$thisType = 'hddm:anyURI'">
              <xsl:attribute name="{@name}">anyURI</xsl:attribute>
            </xsl:when>
            <xsl:when test="$thisType = 'hddm:boolean'">
              <xsl:attribute name="{@name}">boolean</xsl:attribute>
            </xsl:when>
            <xsl:when test="count(ancestor::node()/xs:simpleType[@name=$thisType]) > 0">
              <xsl:apply-templates select="ancestor::node()/xs:simpleType[@name=$thisType]">
                <xsl:with-param name="attName" select="@name"/>
              </xsl:apply-templates>
            </xsl:when>
            <xsl:when test="count(ancestor::node()/xs:simpleType[@name=$localType]) > 0">
              <xsl:apply-templates select="ancestor::node()/xs:simpleType[@name=$localType]">
                <xsl:with-param name="attName" select="@name"/>
              </xsl:apply-templates>
            </xsl:when>
            <xsl:otherwise>
              <xsl:message terminate="yes">
  Format error: attribute "<xsl:value-of select="name()"/>"
  has unsupported type "<xsl:value-of select="$thisType"/>"
  Did you forget the namespace prefix?
              </xsl:message>
            </xsl:otherwise>
          </xsl:choose>
        </xsl:when>
        <xsl:otherwise>
          <xsl:message terminate="yes">
  Format error: attribute "<xsl:value-of select="name()"/>"
  must either be assigned on of the HDDM types or assigned fixed values
          </xsl:message>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:when>
  </xsl:choose>
</xsl:template>

<xsl:template match="xs:simpleType[@name]">
  <xsl:param name="attName"/>
  <xsl:if test="count(./xs:restriction) != 1">
    <xsl:message terminate="yes">
  Format error: attribute type "<xsl:value-of select="@name"/>"
  must derive by restriction from one of the allowed base types
    </xsl:message>
  </xsl:if>
  <xsl:variable name="thisType">
    <xsl:value-of select="./xs:restriction/@base"/>
  </xsl:variable>
  <xsl:variable name="localType">
    <xsl:value-of select="substring-after(./xs:restriction/@base,':')"/>
  </xsl:variable>
  <xsl:choose>
    <xsl:when test="$thisType = 'hddm:int'">
      <xsl:attribute name="{$attName}">int</xsl:attribute>
    </xsl:when>
    <xsl:when test="$thisType = 'hddm:long'">
      <xsl:attribute name="{$attName}">long</xsl:attribute>
    </xsl:when>
    <xsl:when test="$thisType = 'hddm:float'">
      <xsl:attribute name="{$attName}">float</xsl:attribute>
    </xsl:when>
    <xsl:when test="$thisType = 'hddm:double'">
      <xsl:attribute name="{$attName}">double</xsl:attribute>
    </xsl:when>
    <xsl:when test="$thisType = 'hddm:Particle_t'">
      <xsl:attribute name="{$attName}">string</xsl:attribute>
    </xsl:when>
    <xsl:when test="$thisType = 'hddm:string'">
      <xsl:attribute name="{$attName}">string</xsl:attribute>
    </xsl:when>
    <xsl:when test="$thisType = 'hddm:anyURI'">
      <xsl:attribute name="{$attName}">anyURI</xsl:attribute>
    </xsl:when>
    <xsl:when test="$thisType = 'hddm:boolean'">
      <xsl:attribute name="{$attName}">boolean</xsl:attribute>
    </xsl:when>
    <xsl:when test="count(ancestor::node()/xs:simpleType[@name=$localType]) > 0">
      <xsl:apply-templates select="ancestor::node()/xs:simpleType[@name=$localType]">
        <xsl:with-param name="attName">
          <xsl:value-of select="$attName"/>
        </xsl:with-param>
      </xsl:apply-templates>
    </xsl:when>
    <xsl:otherwise>
      <xsl:message terminate="yes">
  Format error: attribute "<xsl:value-of select="attName"/>"
  has unsupported type "<xsl:value-of select="$thisType"/>"
  Did you forget the namespace prefix?
      </xsl:message>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="xs:attributeGroup">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="text()">
</xsl:template>

</xsl:transform>
