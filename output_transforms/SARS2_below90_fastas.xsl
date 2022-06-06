<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" extension-element-prefixes="exsl">
<xsl:output method="text"/>
<xsl:template match="sample">
<xsl:if test="assay[1]/amplicon[1]/breadth &lt; 90">
<xsl:variable name="SAMPLE">
<xsl:choose><xsl:when test="contains(@name, 'WMTS')">
<xsl:value-of select="substring(@name,1,25)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled-')">
<xsl:value-of select="substring(@name,1,26)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled400')">
<xsl:value-of select="substring(@name,1,29)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled1000')">
<xsl:value-of select="substring(@name,1,30)"/>
</xsl:when><xsl:otherwise>
<xsl:value-of select="@name"/>
</xsl:otherwise></xsl:choose>
</xsl:variable>
<exsl:document method="text" href="{/analysis/@run_name}/{$SAMPLE}_inc.fasta">
<xsl:for-each select="assay">
<xsl:text/>><xsl:value-of select="$SAMPLE"/>
<xsl:text>&#xa;</xsl:text>
<xsl:value-of select="./amplicon[1]/gapfilled_consensus_sequence"/>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
<xsl:text/>
</exsl:document>
</xsl:if>
</xsl:template>
</xsl:stylesheet>
