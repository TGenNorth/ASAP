<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" extension-element-prefixes="exsl">
<xsl:output method="text"/>
<xsl:template match="/analysis">   
<xsl:for-each select="sample">
<xsl:variable name="current_sample" select="@name"/>
<xsl:for-each select="assay[not(@type='presence/absence')]">
<xsl:variable name="current_assay" select="@name"/>
<xsl:choose>
<xsl:when test="@type = 'SNP'">
<xsl:for-each select="amplicon//snp/@name[. != 'unknown' and . != 'position of interest']">
<xsl:sort select="."/>
<xsl:value-of select='$current_sample'/><xsl:text>&#x9;</xsl:text>
<xsl:value-of select="$current_assay"/>smor<xsl:value-of select="."/><xsl:text>&#x9;</xsl:text>
<xsl:choose>
<xsl:when test="../significance/@flag"><xsl:value-of select="../significance/@flag"/><xsl:text>&#x9;</xsl:text>-<xsl:text>&#x9;</xsl:text><xsl:value-of select="../../@reads"/></xsl:when>
<xsl:when test="../significance[not(@flag)]">
<xsl:value-of select='format-number(../snp_call/@percent div 100, "0.####")'/><xsl:text>&#x9;</xsl:text>
<xsl:value-of select="../snp_call/@count"/><xsl:text>&#x9;</xsl:text>
<xsl:value-of select="../../@reads"/>
</xsl:when>
<xsl:otherwise>-<xsl:text>&#x9;</xsl:text>-<xsl:text>&#x9;</xsl:text><xsl:value-of select="../../@reads"/></xsl:otherwise>
</xsl:choose>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
</xsl:when>
<xsl:when test="@type = 'ROI'">
<xsl:for-each select="amplicon//region_of_interest//mutation/@name">
<xsl:sort select="."/>
<xsl:value-of select='$current_sample'/><xsl:text>&#x9;</xsl:text>
<xsl:value-of select="$current_assay"/>smor<xsl:value-of select="."/><xsl:text>&#x9;</xsl:text>
<xsl:choose>
<xsl:when test="../../significance/@flag"><xsl:value-of select="../../significance/@flag"/><xsl:text>&#x9;</xsl:text>-<xsl:text>&#x9;</xsl:text><xsl:value-of select="../../../@reads"/></xsl:when>
<xsl:when test="../../significance[not(@flag)]">
<xsl:value-of select='format-number(../@percent div 100, "0.####")'/><xsl:text>&#x9;</xsl:text>
<xsl:value-of select="../@count"/><xsl:text>&#x9;</xsl:text>
<xsl:value-of select="../../../@reads"/>
</xsl:when>
<xsl:otherwise>-<xsl:text>&#x9;</xsl:text>-<xsl:text>&#x9;</xsl:text><xsl:value-of select="../../../@reads"/></xsl:otherwise>
</xsl:choose>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
</xsl:when>
</xsl:choose>
</xsl:for-each>
</xsl:for-each>
<xsl:text/>
</xsl:template>
</xsl:stylesheet>
