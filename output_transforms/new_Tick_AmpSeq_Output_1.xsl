<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:exsl="http://exslt.org/common" xmlns:str="http://exslt.org/strings" extension-element-prefixes="exsl str">
    <xsl:output method="xhtml" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" omit-xml-declaration="yes" encoding="UTF-8" indent="yes"/>

<!-- SNP Summary Report -->
    <xsl:template match="/analysis">
        <html>
        <head>
            <title>Run Summary for: <xsl:value-of select="@run_name"/></title>
        </head>
        <body>
            <center><h1>ASAP Run Summary for: <xsl:value-of select="@run_name"/></h1></center>
	    <br/>

            <table border="1" width="100%">
	        <tr>
	    	<th rowspan = "2">Sample</th>
	    	<xsl:for-each select="sample[1]">
	    	    <xsl:for-each select="assay">
			<xsl:sort select = "name"/>
	    		<xsl:if test="@type = 'SNP' or @type = 'mixed'">
			    <th colspan = "{count(./amplicon/snp[@name!='unknown'])}" nowrap="true"><a href="{/analysis/@run_name}/{@name}.html"><xsl:value-of select="@name"/>:Known SNPs Present</a></th>
			</xsl:if>
			<xsl:if test="@type = 'ROI'">
			  <th colspan = "{count(./amplicon/snp[@name!='unknown'])}" nowrap="true"><a href="{/analysis/@run_name}/{@name}.html"><xsl:value-of select="@name"/></a></th>
			</xsl:if>
	    	    </xsl:for-each>
	    	</xsl:for-each>
		</tr>
		<xsl:for-each select="sample[1]">
		    <tr>
                    <xsl:apply-templates select="."/>
		    <xsl:for-each select="assay">
	    	        <xsl:if test="@type = 'SNP' or @type = 'mixed'">
		    	    <xsl:for-each select="amplicon/snp">
		    	        <xsl:if test="./@name != 'unknown'">
		    	            <th nowrap="true"><xsl:value-of select="./@position"/><xsl:value-of select="./@reference"/>-><xsl:value-of select="./snp_call"/>
		    		    - <xsl:value-of select="significance"/><xsl:if test="significance/@flag">(<xsl:value-of select="significance/@flag"/>)</xsl:if>
		    		    </th>
				</xsl:if>
		    	    </xsl:for-each>
			</xsl:if>
			<xsl:if test="@type = 'ROI'">
			  <th nowrap="true">N/A</th>
			</xsl:if>
	            </xsl:for-each>
		    </tr>
		</xsl:for-each>
		<xsl:for-each select ="sample">
		    <tr>
		    <td nowrap="true"><a href="{/analysis/@run_name}/{./@name}.html"><xsl:value-of select="@name"/></a></td>
		    <xsl:apply-templates select="."/>
		    <xsl:for-each select="assay">
	    	        <xsl:if test="@type = 'SNP' or @type = 'mixed'">
		    	    <xsl:for-each select="amplicon/snp">
		    	        <xsl:if test="./@name != 'unknown'">
		    		    <td align = "center"><xsl:value-of select='format-number(./snp_call/@percent *0.01, "#.###")' />
		    		    </td>
				</xsl:if>
		    	    </xsl:for-each>
			</xsl:if>
			<xsl:if test="@type = 'ROI'">
			  <td align="center">N/A</td>
			</xsl:if>
		    </xsl:for-each>
		    </tr>
		</xsl:for-each>
            </table>
	    <br />
            <a href="{@run_name}_counts.html">Click here for read counts</a>
        </body>
        </html>

<!-- Read Count Summary -->
    <exsl:document method="html" href="{@run_name}_counts.html">
        <html>
        <head>
            <title>Run Summary for: <xsl:value-of select="@run_name"/></title>
        </head>
        <body>
        	<center><h1>ASAP Read Count Summary for: <xsl:value-of select="@run_name"/></h1></center>
	        <br/>
	    	<em>Number of reads aligning to each assay</em>
            <table border="1" cellpadding="3">
	    		<tr>
	    		<th>Sample</th>
	    		<xsl:for-each select="sample[1]">
	    		<xsl:for-each select="assay">
	    		    <th nowrap="true"><xsl:value-of select='@name'/></th>
	    		</xsl:for-each>
	    		</xsl:for-each>
	    		</tr>
                <xsl:for-each select="sample">
                    <tr>
                        <td><a href="{/analysis/@run_name}/{./@name}.html"><xsl:value-of select="@name"/></a></td>
			    		<xsl:for-each select="assay">
			    		    <td align="center">
			    		        <xsl:value-of select="sum(.//amplicon/@reads)"/>
			    		    </td>
			    		</xsl:for-each>
                    </tr>
                    <xsl:apply-templates select="."/>
                </xsl:for-each>
            </table>
	    	<br />
	    	<br />
	    	<a href="{@run_name}.html">Click here for SNP report</a>
        </body>
        </html>
    </exsl:document>
    </xsl:template>

    <xsl:template match="sample">
	<exsl:document method="html" href="{/analysis/@run_name}/{@name}.html">
	    <html>
	    <head>
	    	<title>ASAP Results for Sample: <xsl:value-of select="@name"/></title>
            <xsl:text disable-output-escaping="yes"><![CDATA[<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/2.1.0/Chart.min.js"></script>]]></xsl:text>
	    	<style>
	    	    .ampGraph {
	    	        position: fixed;
	    	        top: 0;
	    	        right: 0;
	    	        bottom: 0;
	    	        left: 0;
	    	        background: rgba(80,80,80,0.8);
	    	        z-index: 99999;
	    	        opacity: 0;
	    	        -webkit-transition: opacity 400ms ease-in;
	                -moz-transition: opacity 400ms ease-in;
	                transition: opacity 400ms ease-in;
	                pointer-events: none;
	    	    }
	    	    .ampGraph:target {
			opacity:1;
			pointer-events: auto;
		    }
		    .ampGraph <xsl:text disable-output-escaping="yes"><![CDATA[>]]></xsl:text> div {
			width: 95vw;
			height: 60vh;
			position: relative;
			margin: 10% auto;
			padding: 5px 20px 13px 20px;
			border-radius: 10px;
			background: #fff;
			background: -moz-linear-gradient(#fff, #999);
			background: -webkit-linear-gradient(#fff, #999);
			background: -o-linear-gradient(#fff, #999);
	            }
		    .close {
			background: #606061;
			color: #FFFFFF;
			line-height: 25px;
			position: absolute;
			right: -12px;
			text-align: center;
			top: -10px;
			width: 24px;
		        text-decoration: none;
			font-weight: bold;
			-webkit-border-radius: 12px;
			-moz-border-radius: 12px;
			border-radius: 12px;
			-moz-box-shadow: 1px 1px 3px #000;
			-webkit-box-shadow: 1px 1px 3px #000;
			box-shadow: 1px 1px 3px #000;
		    }
		    .close:hover { background: #00d9ff; }
		    .ampCanvas {
		        overflow-x: auto;
		    }
	    	</style>
	    </head>
	    <body>
	        <center><h1>ASAP Results for Sample: <xsl:value-of select="@name"/></h1></center>
	        <br />
	        Total reads: <xsl:value-of select="@mapped_reads + @unmapped_reads + @unassigned_reads"/><br/>
	        Mapped reads: <xsl:value-of select="@mapped_reads"/><br/>
	        Unmapped reads: <xsl:value-of select="@unmapped_reads + @unassigned_reads"/><br/>
	        <br/>
	    	<table border="2" rules="rows">
	    	    <tr><th colspan="2">Clinical Summary for Sample: <xsl:value-of select="@name"/></th></tr>
	    	    <xsl:for-each select=".//significance">
                    <xsl:if test="not(./@flag)">
	    		    <tr>
	    		        <td><xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;]]></xsl:text></td>
	    		        <td><xsl:value-of select="."/></td>
	    		    </tr>
                    </xsl:if>
	            </xsl:for-each>
	    	</table>
	    	<br />
	    	<br />
	    	<h3>Presence/absence assays for sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Assay Name</th>
	    		<th># of Reads</th>
	    		<th>Coverage Breadth</th>
	    		<th>Significance</th>
	    		<th>SNPs found(% reads containing SNP)</th>
	    		</tr>
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@type = 'presence/absence'">
	    		    <tr>
	    		        <td><xsl:value-of select="@name"/></td>
	    		        <td><xsl:value-of select="amplicon/@reads"/></td>
	    		        <td><xsl:value-of select='format-number(amplicon/breadth, "##.##")'/>%</td>
	    		        <td><xsl:value-of select="amplicon/significance"/><xsl:if test="amplicon/significance/@flag"> (<xsl:value-of select="amplicon/significance/@flag"/>)</xsl:if></td>
	    		        <td><xsl:for-each select="amplicon/snp">
	    		            <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)<br/>
	    		        </xsl:for-each></td>
	    		    </tr>
	    		    </xsl:if>
	    		</xsl:for-each>
	    	</table>
	    	<br />
	    	<h3>SNP assays for sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Assay Name</th>
	    		<th># of Reads</th>
	    		<th>Coverage Breadth</th>
	    		<th>Known SNPs(% reads containing SNP) - Significance</th>
	    		<th>Other SNPs found(% reads containing SNP)</th>
	    		</tr>
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@type = 'SNP' or @type = 'mixed'">
	    		    <tr>
	    		        <td><xsl:value-of select="@name"/></td>
	    		        <td><xsl:value-of select="amplicon/@reads"/></td>
	    		        <xsl:if test="amplicon/@reads &gt; 0">
		    		        <td><xsl:value-of select='format-number(amplicon/breadth, "##.##")'/>%</td>
		    		        <td><xsl:for-each select="amplicon/snp">
		    		            <xsl:if test="./@name != 'unknown'">
		    		                <xsl:value-of select="./@position"/><xsl:value-of select="./@reference"/>-><xsl:value-of select="./snp_call"/>(<xsl:value-of select='format-number(./snp_call/@percent, "##.##")'/>%)
		    		                - <xsl:value-of select="significance"/><xsl:if test="significance/@flag">(<xsl:value-of select="significance/@flag"/>)</xsl:if>
		    		                <br/>
		    		            </xsl:if>
		    		        </xsl:for-each></td>
		    		        <td><xsl:for-each select="amplicon/snp">
		    		            <xsl:if test="./@name = 'unknown'">
		    		                <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)<br/>
		    		            </xsl:if>
		    		        </xsl:for-each></td>
	    		        </xsl:if>
	    		    </tr>
	    		    </xsl:if>
	    		</xsl:for-each>
	    	</table>
	    	<br />
	    	
                <h3>Gene Variant assays for sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Assay Name</th>
	    		<th>Gene variant - # of Reads - Coverage breadth - Significance</th>
	    		</tr>
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@type = 'gene variant'">
	    		    <tr>
	    		        <td><xsl:value-of select="@name"/></td>
		    		    <td><xsl:for-each select="amplicon">
                                    <xsl:sort select="@reads" data-type="number" order="descending"/>
	    		            <xsl:if test="@reads &gt; 0">
	    		                <xsl:value-of select="@variant"/> - <strong><xsl:value-of select="@reads"/></strong> - <xsl:value-of select='format-number(breadth, "##.##")'/>%
	    		                <xsl:if test="significance">
	    		                    - <xsl:value-of select="significance"/><xsl:if test="significance/@flag"> (<xsl:value-of select="significance/@flag"/>)</xsl:if>
	    		               </xsl:if>
                                <br />
		    		        </xsl:if>
		    		    </xsl:for-each></td>
	    		    </tr>
	    		    </xsl:if>
	    		</xsl:for-each>
	    	</table>
		<!--allele tables-->
		<xsl:for-each select="assay">
		  <xsl:if test="@type='mixed' or @type='ROI'">
		    <h3>Alleles for Region of Interest: <xsl:value-of select="@name"/></h3>
	    	    <table border="2" width="100%">
	    	      <tr>
	    		<th># of Reads</th>
	    		<th>% Frequency</th>
			<th>Allele Sequences</th>
	    	      </tr>
          <xsl:choose>
            <xsl:when test="amplicon/region_of_interest/allele_sequence">
		      <xsl:for-each select="amplicon/region_of_interest/allele_sequence">
          			<tr>
          			  <td><xsl:value-of select="@count"/></td>
          			  <td><xsl:value-of select='format-number(@percent, "##.##")'/>%</td>
          			  <td><a href="{@hash}.html"><xsl:value-of select="."/></a></td>
          			</tr>
		      </xsl:for-each>
        </xsl:when>
        <xsl:otherwise><tr><td>Not Found</td><td></td><td></td></tr></xsl:otherwise>
      </xsl:choose>
		    </table>
		  </xsl:if>
		</xsl:for-each>
	      </body>
	    </html>
	</exsl:document>
    </xsl:template>

</xsl:stylesheet>
