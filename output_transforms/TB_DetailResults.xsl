<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" extension-element-prefixes="exsl">
    <xsl:output method="xhtml" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" omit-xml-declaration="yes" encoding="UTF-8" indent="yes"/>
    <xsl:template match="/analysis">
        <html>
        <head>
            <title>Run Summary for: <xsl:value-of select="@run_name"/></title>
        </head>
        <body>
        	<center><h1>TB Detailed ASAP Run Summary for: <xsl:value-of select="@run_name"/></h1></center>
	        <br/>
	    	<em>Number of reads aligning to each assay</em>
            <table border="1" width="100%">
	    		<tr>
	    		<th>Sample</th>
	    		<xsl:for-each select="sample[1]">
	    		<xsl:for-each select="assay">
	    		    <th><xsl:value-of select='@name'/></th>
	    		</xsl:for-each>
	    		</xsl:for-each>
	    		</tr>
                <xsl:for-each select="sample">
                    <tr>
                        <td><a href="{/analysis/@run_name}/{./@name}_details.html"><xsl:value-of select="@name"/></a></td>
			    		<xsl:for-each select="assay">
			    		    <td align="center"><xsl:for-each select="amplicon">
			    		        <xsl:value-of select="./@reads"/>
			    		    </xsl:for-each></td>
			    		</xsl:for-each>
                    </tr>
                    <xsl:apply-templates select="."/>
                </xsl:for-each>
            </table>
        </body>
        </html>
	</xsl:template>
	<xsl:template match="sample">
	<exsl:document method="html" href="{/analysis/@run_name}/{@name}_details.html">
	    <html>
	    <head>
	    	<title>Detailed Results for Sample: <xsl:value-of select="@name"/></title>
	    </head>
	    <body>
	        <center><h1>Detailed ASAP Results for Sample: <xsl:value-of select="@name"/></h1></center>
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
	    	<h3><em>M. tuberculosis</em> identification for sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Assay Name</th>
	    		<th># of Reads</th>
	    		<th>Coverage Breadth</th>
	    		<th>Significance</th>
	    		<th>SNPs found(% reads containing SNP)</th>
	    		</tr>
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@name = 'IS6110'">
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
	    	<h3>Drug resistance assays for sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Assay Name</th>
	    		<th># of Reads</th>
	    		<th>Known Mutations(% reads containing mutation) - Significance</th>
	    		<th>Other SNPs found(% reads containing SNP)</th>
	    		</tr>
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@type = 'SNP' or @type = 'ROI' or @type = 'mixed'">
	    		    <tr>
	    		        <td><xsl:value-of select="@name"/></td>
	    		        <td><xsl:value-of select="amplicon/@reads"/></td>
	    		        <xsl:if test="amplicon/@reads &gt; 0">
		    		        <td>
		    		        <xsl:for-each select="amplicon/snp">
		    		            <xsl:if test="./@name != 'unknown'">
		    		                <xsl:value-of select="./@name"/>(<xsl:value-of select='format-number(./snp_call/@percent, "##.##")'/>%)
		    		                <xsl:if test="significance"> - <xsl:value-of select="significance"/><xsl:if test="significance/@flag">(<xsl:value-of select="significance/@flag"/>)</xsl:if></xsl:if>
		    		                <br/>
		    		            </xsl:if>
		    		        </xsl:for-each>
		    		        <xsl:for-each select="amplicon/region_of_interest">
		    		            <xsl:for-each select="mutation">
		    		                <xsl:value-of select="./@name"/>(<xsl:value-of select='format-number(./@percent, "##.##")'/>%)
	    		                    <xsl:if test="../significance and ./@percent &gt; 10"> - <xsl:value-of select="../significance"/><xsl:if test="../significance/@flag">(<xsl:value-of select="../significance/@flag"/>)</xsl:if></xsl:if>
		    		                <br/>
		    		            </xsl:for-each>
		    		        </xsl:for-each>
		    		        </td>
		    		        <td>
		    		        <xsl:for-each select="amplicon/snp">
		    		            <xsl:if test="./@name = 'unknown'">
		    		                <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)<br/>
		    		            </xsl:if>
		    		        </xsl:for-each>
		    		        </td>
	    		        </xsl:if>
	    		    </tr>
	    		    </xsl:if>
	    		</xsl:for-each>
	    	</table>
	    	<br />
	    </body>
	    </html>
	</exsl:document>
    </xsl:template>
</xsl:stylesheet>
