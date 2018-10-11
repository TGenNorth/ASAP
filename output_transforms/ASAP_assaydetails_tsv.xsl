<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" extension-element-prefixes="exsl">
<xsl:output method="text"/>

<xsl:template name="compute_proportion">
  <xsl:param name="value"/>
  <xsl:param name="total"/>
  <xsl:if test="$value">
      <xsl:value-of select="format-number($value div $total, '#.###')" />
    </xsl:if>
    <xsl:if test="not($value)">
      <xsl:text>0</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="/analysis">

<xsl:for-each select="sample[1]/assay[@type='SNP']">
<xsl:variable name="ASSAY_NAME" select="@name"/>
<exsl:document method="text" href="{/analysis/@run_name}/{$ASSAY_NAME}.tsv">
<xsl:text/>Sample Name<xsl:text>&#x9;</xsl:text>#Reads<xsl:text>&#x9;</xsl:text>Coverage Breadth<xsl:text>&#x9;</xsl:text>SNP Position<xsl:text>&#x9;</xsl:text>Depth<xsl:text>&#x9;</xsl:text>Reference<xsl:text>&#x9;</xsl:text>A<xsl:text>&#x9;</xsl:text>C<xsl:text>&#x9;</xsl:text>G<xsl:text>&#x9;</xsl:text>T<xsl:text>&#x9;</xsl:text>_<xsl:text>&#x9;</xsl:text>Significance<xsl:text>&#xa;</xsl:text>
<xsl:for-each select="//sample/assay[@name = $ASSAY_NAME]">
<xsl:variable name="SAMPLE_NAME" select="../@name"/>
<xsl:variable name="NUM_READS" select="amplicon/@reads"/>
<xsl:variable name="BREADTH" select="format-number(amplicon/breadth, '##.##')"/>
<xsl:for-each select="amplicon/snp[@name != 'unknown' and not(@position = preceding-sibling::snp/@position)]">
<xsl:variable name="DEPTH" select="@depth"/>
<xsl:value-of select="$SAMPLE_NAME"/><xsl:text>&#x9;</xsl:text><xsl:value-of select="$NUM_READS"/><xsl:text>&#x9;</xsl:text><xsl:value-of select="$BREADTH"/>%<xsl:text>&#x9;</xsl:text><xsl:value-of select="@position"/><xsl:text>&#x9;</xsl:text><xsl:value-of select="$DEPTH"/><xsl:text>&#x9;</xsl:text><xsl:value-of select="@reference"/><xsl:text>&#x9;</xsl:text><xsl:call-template name="compute_proportion"><xsl:with-param name="value" select="base_distribution/@A"/><xsl:with-param name="total" select="$DEPTH"/></xsl:call-template><xsl:text>&#x9;</xsl:text><xsl:call-template name="compute_proportion"><xsl:with-param name="value" select="base_distribution/@C"/><xsl:with-param name="total" select="$DEPTH"/></xsl:call-template><xsl:text>&#x9;</xsl:text><xsl:call-template name="compute_proportion"><xsl:with-param name="value" select="base_distribution/@G"/><xsl:with-param name="total" select="$DEPTH"/></xsl:call-template><xsl:text>&#x9;</xsl:text><xsl:call-template name="compute_proportion"><xsl:with-param name="value" select="base_distribution/@T"/><xsl:with-param name="total" select="$DEPTH"/></xsl:call-template><xsl:text>&#x9;</xsl:text><xsl:call-template name="compute_proportion"><xsl:with-param name="value" select="base_distribution/@_"/><xsl:with-param name="total" select="$DEPTH"/></xsl:call-template><xsl:text>&#x9;</xsl:text><xsl:if test="not(significance/@flag)"><xsl:value-of select="significance"/></xsl:if><xsl:text>&#xa;</xsl:text>
</xsl:for-each>
</xsl:for-each>
</exsl:document>
</xsl:for-each>

<xsl:for-each select="sample[1]/assay[@type='presence/absence']">
<xsl:variable name="ASSAY_NAME" select="@name"/>
<exsl:document method="text" href="{/analysis/@run_name}/{$ASSAY_NAME}.tsv">
<xsl:text/>Sample Name<xsl:text>&#x9;</xsl:text>#Reads<xsl:text>&#x9;</xsl:text>Coverage Breadth<xsl:text>&#x9;</xsl:text>Significance<xsl:text>&#xa;</xsl:text>
<xsl:for-each select="//sample/assay[@name = $ASSAY_NAME]">
<xsl:variable name="SAMPLE_NAME" select="../@name"/>
<xsl:variable name="NUM_READS" select="amplicon/@reads"/>
<xsl:variable name="BREADTH" select="format-number(amplicon/breadth, '##.##')"/>
<xsl:value-of select="$SAMPLE_NAME"/><xsl:text>&#x9;</xsl:text><xsl:value-of select="$NUM_READS"/><xsl:text>&#x9;</xsl:text><xsl:value-of select="$BREADTH"/>%<xsl:text>&#x9;</xsl:text><xsl:if test="not(signifcance/@flag)"><xsl:value-of select="significance"/></xsl:if><xsl:text>&#xa;</xsl:text>
</xsl:for-each>
</exsl:document>
</xsl:for-each>

</xsl:template>
</xsl:stylesheet>
