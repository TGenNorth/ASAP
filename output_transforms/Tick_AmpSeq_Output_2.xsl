<?xml version="1.0" encoding="UTF-8"?>

<xsl:stylesheet
  version="2.0"
  xmlns:asap="http://pathogen.tgen.org/ASAP/functions"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:exsl="http://exslt.org/common"
  extension-element-prefixes="asap exsl"
>

  <xsl:output
      method="xhtml"
      doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN"
      doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd"
      omit-xml-declaration="yes"
      encoding="UTF-8"
      indent="yes"
  />


<xsl:template match="/">
  <xsl:for-each select="analysis/assay[@type='SNP' or @type='mixed']">
    <xsl:variable name="currentAssay" select="."/>
    <exsl:document method="html" href="{../@run_name}/{@name}.html">
    <html>

      <!-- set title -->
      <h1>ASAP Assay Summary: <xsl:value-of select="@name"/></h1>

    <style>
    table, th, td {
     border: 1px solid black;
     border-collapse: collapse;
      }
   </style>
   <body>
      <table style = "width:100%, border: 1px solid black">


<!-- _______SET COLUMN HEADERS ________________________________ -->
        <tr>
          <th>Samples</th>

          <!-- For each snp name print it as a column header -->
          <xsl:variable name="colHeaders" select="asap:distinct-values(sample/amplicon/snp/@name)"/>
          <xsl:for-each select="$colHeaders">
           <xsl:sort select="."/>
           <th><xsl:value-of select="."/></th>
        </xsl:for-each>
        </tr>
  <!-- _____________________________________________________________ -->


  <!-- get the column name list... again. This is done to give the complete list. Before this we were getting a list
 of only columns the sample DID have. Not what they do AND do not have. -->
  <xsl:variable name="colHeaders" select="asap:distinct-values($currentAssay/sample/amplicon/snp/@name)"/>

  <!-- For each sample -->
        <xsl:for-each select="$currentAssay/sample">
          <tr>
            <!-- print the sample name -->
            <td><a href="{@name}.html"><xsl:value-of select="@name"/></a></td>
            <!-- store name in currentSamp for later -->
            <xsl:variable name="currentSamp" select="@name"/>
            <!-- for every snp name in the sample -->
            <xsl:for-each select="$colHeaders">
              <xsl:sort select="."/>
              <!-- save the snp name -->
              <xsl:variable name="currentSnp" select="."/>
              <!-- check to see if the current sample has the current snp (count that is no 0)
              Warning: there should only be 1 snp that matches this. If thats not true this
              will be wrong. -->
              <xsl:variable name="gotValue" select="count($currentAssay/sample[@name = $currentSamp]/amplicon/snp[@name = $currentSnp])"/>

              <!-- choose is an if else statement -->
              <xsl:choose>
                <!-- when we have  a snp in this sample that matches the current one we are looking for.-->
                <xsl:when test="$gotValue > 0">
                      <!-- print the snp's snip_call value as a percent-->
                    <td><xsl:value-of select="format-number($currentAssay/sample[@name = $currentSamp]/amplicon/snp[@name = $currentSnp]/snp_call/@percent div 100,'0.####')"/></td>
                </xsl:when>
                  <!-- otherwise print 0 -->
                <xsl:otherwise>
                  <td>0.0000</td>
                </xsl:otherwise>
              </xsl:choose>
            </xsl:for-each>
          </tr>
        </xsl:for-each>

      </table>
    </body>
  </html>
</exsl:document>
</xsl:for-each>
</xsl:template>
</xsl:stylesheet>
