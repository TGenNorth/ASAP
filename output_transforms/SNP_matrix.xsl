<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:asap="http://pathogen.tgen.org/ASAP/functions" xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:exsl="http://exslt.org/common" xmlns:str="http://exslt.org/strings" extension-element-prefixes="exsl str asap">
<xsl:output method="text" />
<xsl:template match="/analysis">
<xsl:text>snp&#x9;mutation</xsl:text>
<xsl:for-each select="sample">
<xsl:variable name="SAMPLE">
<xsl:choose><xsl:when test="contains(@name, 'WMTS')">
<xsl:value-of select="substring(@name,18,8)"/>
</xsl:when><xsl:when test="contains(@name, 'ARTIC')">
<xsl:value-of select="substring(@name,19,8)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled-')">
<xsl:value-of select="substring(@name,19,9)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled400')">
<xsl:value-of select="substring(@name,22,8)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled1000')">
<xsl:value-of select="substring(@name,23,8)"/>
</xsl:when><xsl:otherwise>
<xsl:value-of select="@name"/>
</xsl:otherwise></xsl:choose>
</xsl:variable>
<xsl:text>&#x9;</xsl:text><xsl:value-of select="$SAMPLE"/>
</xsl:for-each>
<xsl:text>&#xa;</xsl:text>
<xsl:for-each select="asap:distinct-values(//snp/@position)">
<xsl:sort data-type="number" select="."/>
<xsl:variable name="POSITION" select="."/>
<xsl:value-of select="//snp[@position=$POSITION]/@reference"/><xsl:value-of select="$POSITION"/><xsl:value-of select="//snp[@position=$POSITION]/snp_call"/>
<xsl:text>&#x9;</xsl:text><xsl:value-of select="//snp[@position=$POSITION]/@name"/>
<xsl:for-each select="//sample">
<xsl:choose>
<xsl:when test=".//snp/@position = $POSITION">
<xsl:variable name="SNP" select=".//snp[@position=$POSITION]"/>
<xsl:text>&#x9;</xsl:text><xsl:value-of select="$SNP/snp_call/@count"/>/<xsl:value-of select="$SNP/@depth"/>
</xsl:when>
<xsl:otherwise>
<xsl:text>&#x9;0.0</xsl:text>
</xsl:otherwise>
</xsl:choose>
</xsl:for-each>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
<xsl:text/>
</xsl:template>
</xsl:stylesheet>
