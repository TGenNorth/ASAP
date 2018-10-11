<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" extension-element-prefixes="exsl">
<xsl:output method="text"/>
<xsl:template match="/analysis/assay">    
<exsl:document method="text" href="{@name}.csv">
<xsl:for-each select="sample">
<xsl:if test="amplicon[@reads&gt;0]">
<xsl:text/><xsl:value-of select="@name"/>,<xsl:value-of select="amplicon/depths"/>
<xsl:text>&#xa;</xsl:text>
<xsl:text>&#xa;</xsl:text>
<xsl:text/><xsl:value-of select="@name"/>,<xsl:value-of select="amplicon/discards"/>
<xsl:text>&#xa;</xsl:text>
<xsl:text>&#xa;</xsl:text>
<xsl:text>&#xa;</xsl:text>
</xsl:if>
</xsl:for-each>
<xsl:text/>
</exsl:document>
</xsl:template>
</xsl:stylesheet>
