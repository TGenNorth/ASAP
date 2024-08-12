<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:exsl="http://exslt.org/common" xmlns:str="http://exslt.org/strings" extension-element-prefixes="exsl str">
<xsl:output method="text" />
<xsl:text/>
<xsl:template match="/analysis">
<xsl:for-each select="sample">
<xsl:if test="assay[1]/amplicon[1]/breadth &gt; 75 and assay[1]/amplicon[1]/average_depth &gt; 30">
<xsl:variable name="SAMPLE">
<xsl:choose><xsl:when test="contains(@name, 'WMTS')">
<xsl:value-of select="substring(@name,18,8)"/>
</xsl:when><xsl:when test="contains(@name, 'ARTIC')">
<xsl:value-of select="substring(@name,19,8)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled-TG1')">
<xsl:value-of select="substring(@name,19,9)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled-')">
<xsl:value-of select="substring(@name,19,8)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled400')">
<xsl:value-of select="substring(@name,22,8)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled1000')">
<xsl:value-of select="substring(@name,23,8)"/>
</xsl:when><xsl:otherwise>
<xsl:value-of select="@name"/>
</xsl:otherwise></xsl:choose>
</xsl:variable>
<xsl:for-each select="assay[1]/amplicon[1]/snp">
<xsl:if test="@depth &gt;= 10 and snp_call/@percent &gt;= 80">
<xsl:value-of select="$SAMPLE"/><xsl:text>&#x9;</xsl:text><xsl:value-of select="@position"/><xsl:text>&#x9;</xsl:text><xsl:value-of select="@reference"/><xsl:text>&#x9;</xsl:text><xsl:value-of select="snp_call"/><xsl:text>&#xa;</xsl:text>
</xsl:if>
</xsl:for-each>
</xsl:if>
</xsl:for-each>
<xsl:text/>
</xsl:template>
</xsl:stylesheet>
