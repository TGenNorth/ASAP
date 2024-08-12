<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" extension-element-prefixes="exsl">
<xsl:output method="text"/>
<xsl:template match="/analysis">    
<xsl:text/>sample_name,assay,codon_name,reference,call,distribution<xsl:text>&#xa;</xsl:text>
<xsl:for-each select="sample">
<xsl:variable name="SAMPLE" select="@name"/>
<xsl:for-each select="assay">
<xsl:variable name="ASSAY" select="@name"/>
<xsl:for-each select=".//region_of_interest">
<xsl:text/><xsl:value-of select="$SAMPLE"/>,<xsl:value-of select="$ASSAY"/>,<xsl:value-of select="@name"/>,<xsl:value-of select="@aa_reference"/>,<xsl:value-of select="amino_acid_sequence"/>,<xsl:for-each select="aa_sequence_distribution/@*"><xsl:value-of select="name()"/>=<xsl:value-of select="."/><xsl:if test="position() != last()">|</xsl:if></xsl:for-each>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
</xsl:for-each>
</xsl:for-each>
<xsl:text/>
</xsl:template>
</xsl:stylesheet>
