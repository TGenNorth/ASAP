<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" extension-element-prefixes="exsl">
<xsl:output method="text"/>
<xsl:template match="/">
  <xsl:for-each select="analysis/assay_counts/assay">
    <xsl:variable name="currentAssay" select="."/>
    <exsl:document method="text" href="./fasta/{$currentAssay/@name}.fasta" >
      <xsl:for-each select="allele/sample">
        <xsl:sort select="@count" order="descending" data-type="number"/>
	       <xsl:text/>><xsl:value-of select="@name"/>_<xsl:value-of select="@count"/>
	        <xsl:text>&#xa;</xsl:text>
	         <xsl:value-of select="../@sequence"/>
	          <xsl:text>&#xa;</xsl:text>
      </xsl:for-each>
       <xsl:text>&#xa;</xsl:text>
    </exsl:document>
  </xsl:for-each>
</xsl:template>
</xsl:stylesheet>
