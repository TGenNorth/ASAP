<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" extension-element-prefixes="exsl">
<xsl:output method="text"/>
<xsl:template match="/analysis">   
<xsl:text/>Target<xsl:for-each select="sample"><xsl:text>&#x9;</xsl:text><xsl:value-of select="@name"/></xsl:for-each>
<xsl:text>&#xa;</xsl:text>
<xsl:for-each select="sample[1]">
<xsl:for-each select="assay">
<xsl:value-of select='@name'/>
<xsl:variable name="CURRENT" select="@name"/>
<xsl:for-each select="//sample">
<xsl:text>&#x9;</xsl:text><xsl:value-of select="./assay[@name=$CURRENT]/amplicon/@reads"/>
</xsl:for-each>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
</xsl:for-each>
<xsl:text/>
</xsl:template>
</xsl:stylesheet>
