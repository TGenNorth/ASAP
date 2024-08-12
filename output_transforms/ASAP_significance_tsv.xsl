<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" extension-element-prefixes="exsl">
<xsl:output method="text"/>
<xsl:template match="/analysis">   
<xsl:for-each select="sample">
<xsl:variable name="current_sample" select="@name"/>
<xsl:for-each select="assay">
<xsl:variable name="current_assay" select="@name"/>
<xsl:for-each select="descendant::significance[not(@flag)]">
<xsl:value-of select='$current_sample'/><xsl:text>&#x9;</xsl:text>
<xsl:value-of select="$current_assay"/><xsl:text>&#x9;</xsl:text>
<xsl:value-of select="parent::node()/@name"/><xsl:text>&#x9;</xsl:text>
<xsl:value-of select="."/><xsl:text>&#x9;</xsl:text>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
</xsl:for-each>
</xsl:for-each>
<xsl:text/>
</xsl:template>
</xsl:stylesheet>
