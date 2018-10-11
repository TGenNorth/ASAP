<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" extension-element-prefixes="exsl">
<xsl:output method="text"/>
<xsl:template match="/analysis">   
<xsl:text/>Sample<xsl:for-each select="sample[1]"><xsl:for-each select="assay[not(@type='presence/absence')]"><xsl:choose><xsl:when test="@type = 'SNP'"><xsl:for-each select="amplicon//snp/@name[. != 'unknown' and . != 'position of interest']"><xsl:sort select="."/><xsl:text>&#x9;</xsl:text><xsl:value-of select="ancestor::assay/@name"/>/<xsl:value-of select="."/></xsl:for-each></xsl:when><xsl:when test="@type = 'ROI'"><xsl:for-each select="amplicon//region_of_interest//mutation/@name"><xsl:sort select="."/><xsl:text>&#x9;</xsl:text><xsl:value-of select="ancestor::assay/@name"/>/<xsl:value-of select="."/></xsl:for-each></xsl:when></xsl:choose></xsl:for-each></xsl:for-each>
<xsl:text>&#xa;</xsl:text>
<xsl:for-each select="sample">
<xsl:value-of select="@name"/><xsl:for-each select="assay[not(@type='presence/absence')]"><xsl:choose><xsl:when test="@type = 'SNP'"><xsl:for-each select="amplicon//snp/@name[. != 'unknown' and . != 'position of interest']"><xsl:sort select="."/>
<xsl:text>&#x9;</xsl:text><xsl:value-of select="../@depth"/></xsl:for-each></xsl:when><xsl:when test="@type = 'ROI'"><xsl:for-each select="amplicon//region_of_interest//mutation/@name"><xsl:sort select="."/><xsl:text>&#x9;</xsl:text><xsl:value-of select="../../@depth"/></xsl:for-each></xsl:when></xsl:choose></xsl:for-each>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
<xsl:text/>
</xsl:template>
</xsl:stylesheet>
