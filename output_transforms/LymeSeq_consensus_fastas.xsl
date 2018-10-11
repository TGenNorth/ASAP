<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" extension-element-prefixes="exsl">
<xsl:output method="text"/>
<xsl:template match="/analysis//sample">    
<xsl:for-each select="assay">
<xsl:if test=".//amplicon[@reads&gt;0]">
<exsl:document method="text" href="consensus/{@name}_{parent::sample/@name}.fasta">
<xsl:for-each select="amplicon[@reads&gt;0]">
<xsl:choose>
<xsl:when test="@variant"><xsl:text/>><xsl:value-of select="ancestor::sample/@name"/>_<xsl:value-of select="@variant"/></xsl:when>
<xsl:otherwise><xsl:text/>><xsl:value-of select="ancestor::sample/@name"/></xsl:otherwise>
</xsl:choose>
<xsl:text>&#xa;</xsl:text>
<xsl:value-of select="consensus_sequence"/>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
<xsl:text/>
</exsl:document>
</xsl:if>
</xsl:for-each>
</xsl:template>
</xsl:stylesheet>
