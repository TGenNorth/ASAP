<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" extension-element-prefixes="exsl">
<xsl:output method="text"/>
<xsl:template match="/analysis">    
<xsl:for-each select="sample">
<xsl:text/>><xsl:value-of select="@name"/>
<xsl:text>&#xa;</xsl:text>
<xsl:value-of select="./assay/amplicon/consensus_sequence"/>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
<xsl:text/>
</xsl:template>
</xsl:stylesheet>
