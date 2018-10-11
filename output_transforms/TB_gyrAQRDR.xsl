<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:exsl="http://exslt.org/common" xmlns:qrdr="http://tb_qrdr.data" extension-element-prefixes="exsl">
    <xsl:output method="html" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" omit-xml-declaration="yes" encoding="UTF-8" indent="yes"/>
    
    <xsl:key name="qrdr-lookup" match="qrdr:variant" use="qrdr:sequence"/>

    <qrdr:variants>
        <qrdr:variant><qrdr:name>Wildtype</qrdr:name><qrdr:sequence>GDASIYDS</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>G88A</qrdr:name><qrdr:sequence>ADASIYDS</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>G88A/S95T</qrdr:name><qrdr:sequence>ADASIYDT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>G88C/S95T</qrdr:name><qrdr:sequence>CDASIYDT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>G88C/D94G/S95T</qrdr:name><qrdr:sequence>CDASIYGT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>G88C/D94N/S95T</qrdr:name><qrdr:sequence>CDASIYNT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>G88V/S95T</qrdr:name><qrdr:sequence>VDASIYDT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D89G</qrdr:name><qrdr:sequence>GGASIYDS</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D89G/S95T</qrdr:name><qrdr:sequence>GGASIYDT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D89N/S95T</qrdr:name><qrdr:sequence>GNASIYDT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D89N/D94A/S95T</qrdr:name><qrdr:sequence>GNASIYAT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D89Y/A90V/S95T</qrdr:name><qrdr:sequence>GYVSIYDT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>A90T/S95T</qrdr:name><qrdr:sequence>GDTSIYDT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>A90V</qrdr:name><qrdr:sequence>GDVSIYDS</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>A90V/S95T</qrdr:name><qrdr:sequence>GDVSIYDT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>A90V/S91P/S95T</qrdr:name><qrdr:sequence>GDVPIYDT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>A90V/D94A/S95T</qrdr:name><qrdr:sequence>GDVSIYAT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>A90V/D94G</qrdr:name><qrdr:sequence>GDVSIYGS</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>A90V/D94G/S95T</qrdr:name><qrdr:sequence>GDVSIYGT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>A90V/D94N/S95T</qrdr:name><qrdr:sequence>GDVSIYNT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>A90V/D94Y/S95T</qrdr:name><qrdr:sequence>GDVSIYYT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>S91P/S95T</qrdr:name><qrdr:sequence>GDAPIYDT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>S91P/D94A/S95T</qrdr:name><qrdr:sequence>GDAPIYAT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>S91P/D94G/S95T</qrdr:name><qrdr:sequence>GDAPIYGT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>S91P/D94H/S95T</qrdr:name><qrdr:sequence>GDAPIYHT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D94A</qrdr:name><qrdr:sequence>GDASIYAS</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D94A/S95T</qrdr:name><qrdr:sequence>GDASIYAT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D94C/S95T</qrdr:name><qrdr:sequence>GDASIYCT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D94G</qrdr:name><qrdr:sequence>GDASIYGS</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D94G/S95A</qrdr:name><qrdr:sequence>GDASIYGA</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D94G/S95T</qrdr:name><qrdr:sequence>GDASIYGT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D94H/S95T</qrdr:name><qrdr:sequence>GDASIYHT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D94N</qrdr:name><qrdr:sequence>GDASIYNS</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D94N/S95T</qrdr:name><qrdr:sequence>GDASIYNT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D94S/S95T</qrdr:name><qrdr:sequence>GDASIYST</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D94T/S95T</qrdr:name><qrdr:sequence>GDASIYTT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D94V/S95T</qrdr:name><qrdr:sequence>GDASIYVT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>D94Y/S95T</qrdr:name><qrdr:sequence>GDASIYYT</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>S95A</qrdr:name><qrdr:sequence>GDASIYDA</qrdr:sequence></qrdr:variant>
        <qrdr:variant><qrdr:name>S95T</qrdr:name><qrdr:sequence>GDASIYDT</qrdr:sequence></qrdr:variant>
    </qrdr:variants>

    <xsl:template match="qrdr:variants">
        <xsl:param name="sequence"/>
        <xsl:value-of select="key('qrdr-lookup', $sequence)/qrdr:name"/>
    </xsl:template>

    <xsl:variable name="qrdr-top" select="document('')/*/qrdr:variants"/>

    <xsl:template match="/analysis">
        <html>
        <head>
            <title>Run Summary for: <xsl:value-of select="@run_name"/></title>
            <style type="text/css">
                                .table-column-locked {
                                  border-collapse:separate;
                                  border-top: 1px solid black;
                                  border-spacing: 0px;
                                }
                                .div-table-column-locked {
                                  overflow-x: scroll;
                                  margin-left: 11em;
                                  overflow-y: visible;
                                  padding-bottom: 1px;
                                }
                                .table-column-locked td, .table-column-locked th {
                                  margin: 0;
                                  border: 1px solid black;
                                  border-top-width: 0px;
                                  white-space: nowrap;
                                }
                                .table-column-locked th.headcol {
                                  position: absolute;
                                  width: 11em;
                                  left: 0;
                                  top: auto;
                                  border-right: 3px solid black;
                                  border-top-width: 1px;
                                  margin-top: -1px;
                                }

            </style>
        </head>
        <body>
            <center><h1>ASAP Run Summary for: <xsl:value-of select="@run_name"/></h1></center>
	    <br/>
	    <em>Number of reads (% of total) having each QRDR variant</em>
            <div class="div-table-column-locked"><table class="table-column-locked">
	        <tr>
	    	<th class="headcol">Sample</th>
	    	<th>Aligned Reads</th>
	    	<th>Read Depth at QRDR</th>
	    	<xsl:for-each select="document('')/*/qrdr:variants/qrdr:variant">
	    	    <th><xsl:value-of select="qrdr:name"/></th>
	    	</xsl:for-each>
	    	</tr>
                <xsl:for-each select="sample">
                <tr>
                <th class="headcol" nowrap="true"><a href="{/analysis/@run_name}/{./@name}.html"><xsl:value-of select="@name"/></a></th>
		<xsl:for-each select="assay">
                    <xsl:variable name="current_amp"><xsl:copy-of select="./amplicon/*"/></xsl:variable>
		    <td align="center"><xsl:value-of select="amplicon/@reads"/></td>
                    <td><xsl:value-of select="amplicon/region_of_interest/@depth"/></td>
	    	    <xsl:for-each select="document('')/*/qrdr:variants/qrdr:variant">
                        <xsl:variable name="aa_sequence" select="qrdr:sequence"/>
                        <td nowrap="true">
                        <xsl:for-each select="exsl:node-set($current_amp)/region_of_interest/aa_sequence_distribution/@*">
                        <xsl:if test="name()=$aa_sequence and . &gt; ../../@depth * 0.01">
                            <xsl:value-of select="."/> (<xsl:value-of select="format-number(. div ../../@depth * 100, '##.##')"/>%)
                        </xsl:if>
                        </xsl:for-each>
                        </td>
                    </xsl:for-each>
		</xsl:for-each>
                </tr>
                <xsl:apply-templates select="."/>
                </xsl:for-each>
            </table></div>
        </body>
        </html>
    </xsl:template>

    <xsl:template match="sample">
	<exsl:document method="html" href="{/analysis/@run_name}/{@name}.html">
	    <html>
	    <head>
	    	<title>ASAP Results for Sample: <xsl:value-of select="@name"/></title>
	    </head>
	    <body>
	        <center><h1>ASAP Results for Sample: <xsl:value-of select="@name"/></h1></center>
	        <br />
	        Total reads: <xsl:value-of select="@mapped_reads + @unmapped_reads + @unassigned_reads"/><br/>
	        Mapped reads: <xsl:value-of select="@mapped_reads"/><br/>
	        Unmapped reads: <xsl:value-of select="@unmapped_reads + @unassigned_reads"/><br/>
	        <br/>
	    	<br />
	    	<h3>QRDR haplotype analysis for sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th># of Reads</th>
                        <th>Depth at QRDR</th>
	    		<th>Region - Reference - Most common sequence (% reads containing seq) - Significance (# of aa changes)</th>
	    		<th>SNPs found(% reads containing SNP)</th>
	    		</tr>
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@type = 'ROI' or @type = 'mixed' and amplicon/region_of_interest">
	    		    <tr valign="top">
	    		        <td><xsl:value-of select="amplicon/@reads"/></td>
                                <td><xsl:value-of select="amplicon/region_of_interest/@depth"/></td>
	    		        <xsl:if test="amplicon/@reads &gt; 0">
		    		        <td><xsl:for-each select="amplicon/region_of_interest">
		    		            <xsl:value-of select="@region"/> - <xsl:value-of select="@reference"/> - <xsl:value-of select="amino_acid_sequence"/>(<xsl:value-of select='format-number(amino_acid_sequence/@percent, "##.##")'/>%)
	    		                - <xsl:value-of select="significance"/>(<xsl:value-of select="significance/@changes"/>)
		    		            <br/>
		    		            <br/>
                                            <table border="1">
                                                <tr><th>Variant name - AA sequence - # reads (%)</th><th>NT sequence - # reads (%)</th></tr>
                                                <tr valign="top">
                                                <td><xsl:for-each select="aa_sequence_distribution/@*"><xsl:sort select="." data-type="number" order="descending"/>
                                                <xsl:if test=". &gt; ../../@depth * 0.01"><xsl:apply-templates select="$qrdr-top"><xsl:with-param name="sequence" select="name()"/></xsl:apply-templates> - <xsl:value-of select="name()"/> - <xsl:value-of select="."/> (<xsl:value-of select="format-number(. div ../../@depth * 100, '##.##')"/>%)<br /></xsl:if></xsl:for-each></td>
                                                <td><xsl:for-each select="nt_sequence_distribution/@*"><xsl:sort select="." data-type="number" order="descending"/>
                                                <xsl:if test=". &gt; ../../@depth * 0.01"><xsl:value-of select="name()"/> - <xsl:value-of select="."/> (<xsl:value-of select="format-number(. div ../../@depth * 100, '##.##')"/>%)<br /></xsl:if></xsl:for-each></td>
                                                </tr>
                                            </table>
		    		        </xsl:for-each></td>
		    		        <td><xsl:for-each select="amplicon/snp">
		    		            <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)<br/>
		    		        </xsl:for-each></td>
	    		        </xsl:if>
	    		    </tr>
	    		    </xsl:if>
	    		</xsl:for-each>
	    	</table>
	    </body>
	    </html>
	</exsl:document>
    </xsl:template>
</xsl:stylesheet>
