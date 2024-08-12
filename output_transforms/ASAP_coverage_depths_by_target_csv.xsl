<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" extension-element-prefixes="exsl">
<xsl:output method="text"/>
<xsl:template match="/analysis">    
<xsl:for-each select="sample[1]/assay">
<xsl:variable name="ASSAY" select="@name"/>
<xsl:text/><xsl:value-of select="$ASSAY"/>-positions,<xsl:value-of select="amplicon[1]/ref_positions"/>
<xsl:text>&#xa;</xsl:text>
<xsl:for-each select="//assay[@name=$ASSAY]">
<xsl:text/><xsl:value-of select="../@name"/>,<xsl:value-of select="amplicon[1]/depths"/>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
<xsl:text/>
</xsl:template>
</xsl:stylesheet>
