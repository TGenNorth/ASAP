<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:asap="http://pathogen.tgen.org/ASAP/functions" xmlns:exsl="http://exslt.org/common" xmlns:str="http://exslt.org/strings" extension-element-prefixes="exsl str asap">
    <xsl:output method="xhtml" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" omit-xml-declaration="yes" encoding="UTF-8" indent="yes"/>

<!-- Clinical Run Summary -->
    <xsl:template match="/analysis">
        <html>
        <head>
            <title>TB Run Summary for: <xsl:value-of select="@run_name"/></title>
            <style type="text/css">
                .freeze-table {
                    border-spacing: 0;
                    font-size: 14px;
                    padding: 0;
                    border: 1px solid #ccc;
                }
                th {
                    top: 0;
                    position: sticky;
                    background-color: #666;
                    color: #fff;
                    z-index: 20;
                    min-height: 30px;
                    height: 30px;
                    text-align: center;
                    border: 1px solid #fff;
                    white-space: nowrap;
                    padding: 5px;
                }
                .second-header th{
                    top: 30px;
                }
                td {
                    border: 1px solid #ccc;
                    background-color: #fff;
                    white-space: nowrap;
                }
                .fixed-header {
                    z-index: 50;
                }
                .fixed-col-sample {
                    left: 0;
                    position: sticky;
                    width: 240px;
                }
                .fixed-col-breadth {
                    left: 240px;
                    position: sticky;
                    border-right: 3px solid #ccc
                }
            </style>
        </head>
        <body>
        <xsl:variable name="prop_filter" select="sample[1]/@proportion_filter * 100"/>
        <xsl:variable name="mutant_count_filter">
            <xsl:choose>
                <xsl:when test="sample[1]/@mutation_depth_filter"><xsl:value-of select="sample[1]/@mutation_depth_filter"/></xsl:when>
                <xsl:otherwise>0</xsl:otherwise>
            </xsl:choose>
        </xsl:variable>
        	<center><h1>TB ASAP Run Summary for: <xsl:value-of select="@run_name"/></h1></center>
	        <br />
	        <br />
                <!--<div class="div-table-column-locked"><table class="table-column-locked">-->
                <table class="freeze-table">
                    <thead>
	    	    <tr>
	    		<th class="fixed-col-sample fixed-header">Sample</th>
	    		<xsl:for-each select="sample[1]/assay">
	    		    <th nowrap="true"><xsl:value-of select='@name'/></th>
	    		</xsl:for-each>
	    		<xsl:for-each select="sample[1]/assay[starts-with(@name, 'ddn')]/amplicon//snp/@name[. != 'unknown' and . != 'position of interest']">
	    		    <th nowrap="true">
                                ddn <xsl:value-of select="."/>
	    		    </th>
	    		</xsl:for-each>
                        <th>ddn SNPs and INDELs</th>
	    		<xsl:for-each select="sample[1]/assay[starts-with(@name, 'fbiA')]/amplicon//snp/@name[. != 'unknown' and . != 'position of interest']">
	    		    <th nowrap="true">
                                fbiA <xsl:value-of select="."/>
	    		    </th>
	    		</xsl:for-each>
                        <th>fbiA SNPs and INDELs</th>
	    		<xsl:for-each select="sample[1]/assay[starts-with(@name, 'fbiB')]/amplicon//snp/@name[. != 'unknown' and . != 'position of interest']">
	    		    <th nowrap="true">
                                fbiB <xsl:value-of select="."/>
	    		    </th>
	    		</xsl:for-each>
                        <th>fbiB SNPs and INDELs</th>
	    		<xsl:for-each select="sample[1]/assay[starts-with(@name, 'fbiC')]/amplicon//snp/@name[. != 'unknown' and . != 'position of interest']">
	    		    <th nowrap="true">
                                fbiC <xsl:value-of select="."/>
	    		    </th>
	    		</xsl:for-each>
                        <th>fbiC SNPs and INDELs</th>
	    		<xsl:for-each select="sample[1]/assay[starts-with(@name, 'fbiD')]/amplicon//snp/@name[. != 'unknown' and . != 'position of interest']">
	    		    <th nowrap="true">
                                fbiD <xsl:value-of select="."/>
	    		    </th>
	    		</xsl:for-each>
                        <th>fbiD SNPs and INDELs</th>
	    		<xsl:for-each select="sample[1]/assay[starts-with(@name, 'fgd1')]/amplicon//snp/@name[. != 'unknown' and . != 'position of interest']">
	    		    <th nowrap="true">
                                fgd1 <xsl:value-of select="."/>
	    		    </th>
	    		</xsl:for-each>
                        <th>fgd1 SNPs and INDELs</th>
	    		<xsl:for-each select="sample[1]/assay[starts-with(@name, 'fgd2')]/amplicon//snp/@name[. != 'unknown' and . != 'position of interest']">
	    		    <th nowrap="true">
                                fgd2 <xsl:value-of select="."/>
	    		    </th>
	    		</xsl:for-each>
                        <th>fgd2 SNPs and INDELs</th>
	    		</tr>
                        </thead>
                        <xsl:for-each select="sample">
                        <tr>
                        <td class="fixed-col-sample"><a href="{/analysis/@run_name}/{./@name}.html"><xsl:value-of select="@name"/></a></td>
			<xsl:for-each select="assay">
			    <td align="center"><xsl:value-of select="amplicon/@reads"/></td>
                        </xsl:for-each>
                        <xsl:for-each select="assay[starts-with(@name, 'ddn')]/amplicon//snp[@name != 'unknown' and @name != 'position of interest']">
		            <td><xsl:choose>
		            <xsl:when test="../significance/@flag"><font color="lightgray"><em><xsl:value-of select="../significance/@flag"/></em></font></xsl:when>
		            <xsl:when test="./significance/@flag"><font color="lightgray"><xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%) - <em><xsl:value-of select="./significance/@flag"/></em></font></xsl:when>
		            <xsl:when test="./significance[not(@flag)]">
		                <xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)
		            </xsl:when>
		            <xsl:otherwise><!-- SNP not present --><em>-</em></xsl:otherwise>
		            </xsl:choose></td>
                        </xsl:for-each>
	    		<td><xsl:for-each select="assay[starts-with(@name, 'ddn')]/amplicon/snp">
                            <xsl:if test="snp_call/@count &gt;= $mutant_count_filter and snp_call/@percent &gt;= $prop_filter">
	    		    <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/><xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;]]></xsl:text><xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%),<xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;]]></xsl:text>
                            </xsl:if>
	    		</xsl:for-each></td>
                        <xsl:for-each select="assay[starts-with(@name, 'fbiA')]/amplicon//snp[@name != 'unknown' and @name != 'position of interest']">
		            <td><xsl:choose>
		            <xsl:when test="../significance/@flag"><font color="lightgray"><em><xsl:value-of select="../significance/@flag"/></em></font></xsl:when>
		            <xsl:when test="./significance/@flag"><font color="lightgray"><xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%) - <em><xsl:value-of select="./significance/@flag"/></em></font></xsl:when>
		            <xsl:when test="./significance[not(@flag)]">
		                <xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)
		            </xsl:when>
		            <xsl:otherwise><!-- SNP not present --><em>-</em></xsl:otherwise>
		            </xsl:choose></td>
                        </xsl:for-each>
	    		<td><xsl:for-each select="assay[starts-with(@name, 'fbiA')]/amplicon/snp">
                            <xsl:if test="snp_call/@count &gt;= $mutant_count_filter and snp_call/@percent &gt;= $prop_filter">
	    		    <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/><xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;]]></xsl:text><xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%),<xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;]]></xsl:text>
                            </xsl:if>
	    		</xsl:for-each></td>
                        <xsl:for-each select="assay[starts-with(@name, 'fbiB')]/amplicon//snp[@name != 'unknown' and @name != 'position of interest']">
		            <td><xsl:choose>
		            <xsl:when test="../significance/@flag"><font color="lightgray"><em><xsl:value-of select="../significance/@flag"/></em></font></xsl:when>
		            <xsl:when test="./significance/@flag"><font color="lightgray"><xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%) - <em><xsl:value-of select="./significance/@flag"/></em></font></xsl:when>
		            <xsl:when test="./significance[not(@flag)]">
		                <xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)
		            </xsl:when>
		            <xsl:otherwise><!-- SNP not present --><em>-</em></xsl:otherwise>
		            </xsl:choose></td>
                        </xsl:for-each>
	    		<td><xsl:for-each select="assay[starts-with(@name, 'fbiB')]/amplicon/snp">
                            <xsl:if test="snp_call/@count &gt;= $mutant_count_filter and snp_call/@percent &gt;= $prop_filter">
	    		    <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/><xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;]]></xsl:text><xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%),<xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;]]></xsl:text>
                            </xsl:if>
	    		</xsl:for-each></td>
                        <xsl:for-each select="assay[starts-with(@name, 'fbiC')]/amplicon//snp[@name != 'unknown' and @name != 'position of interest']">
		            <td><xsl:choose>
		            <xsl:when test="../significance/@flag"><font color="lightgray"><em><xsl:value-of select="../significance/@flag"/></em></font></xsl:when>
		            <xsl:when test="./significance/@flag"><font color="lightgray"><xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%) - <em><xsl:value-of select="./significance/@flag"/></em></font></xsl:when>
		            <xsl:when test="./significance[not(@flag)]">
		                <xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)
		            </xsl:when>
		            <xsl:otherwise><!-- SNP not present --><em>-</em></xsl:otherwise>
		            </xsl:choose></td>
                        </xsl:for-each>
	    		<td><xsl:for-each select="assay[starts-with(@name, 'fbiC')]/amplicon/snp">
                            <xsl:if test="snp_call/@count &gt;= $mutant_count_filter and snp_call/@percent &gt;= $prop_filter">
	    		    <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/><xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;]]></xsl:text><xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%),<xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;]]></xsl:text>
                            </xsl:if>
	    		</xsl:for-each></td>
                        <xsl:for-each select="assay[starts-with(@name, 'fbiD')]/amplicon//snp[@name != 'unknown' and @name != 'position of interest']">
		            <td><xsl:choose>
		            <xsl:when test="../significance/@flag"><font color="lightgray"><em><xsl:value-of select="../significance/@flag"/></em></font></xsl:when>
		            <xsl:when test="./significance/@flag"><font color="lightgray"><xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%) - <em><xsl:value-of select="./significance/@flag"/></em></font></xsl:when>
		            <xsl:when test="./significance[not(@flag)]">
		                <xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)
		            </xsl:when>
		            <xsl:otherwise><!-- SNP not present --><em>-</em></xsl:otherwise>
		            </xsl:choose></td>
                        </xsl:for-each>
	    		<td><xsl:for-each select="assay[starts-with(@name, 'fbiD')]/amplicon/snp">
                            <xsl:if test="snp_call/@count &gt;= $mutant_count_filter and snp_call/@percent &gt;= $prop_filter">
	    		    <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/><xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;]]></xsl:text><xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%),<xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;]]></xsl:text>
                            </xsl:if>
	    		</xsl:for-each></td>
                        <xsl:for-each select="assay[starts-with(@name, 'fgd1')]/amplicon//snp[@name != 'unknown' and @name != 'position of interest']">
		            <td><xsl:choose>
		            <xsl:when test="../significance/@flag"><font color="lightgray"><em><xsl:value-of select="../significance/@flag"/></em></font></xsl:when>
		            <xsl:when test="./significance/@flag"><font color="lightgray"><xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%) - <em><xsl:value-of select="./significance/@flag"/></em></font></xsl:when>
		            <xsl:when test="./significance[not(@flag)]">
		                <xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)
		            </xsl:when>
		            <xsl:otherwise><!-- SNP not present --><em>-</em></xsl:otherwise>
		            </xsl:choose></td>
                        </xsl:for-each>
	    		<td><xsl:for-each select="assay[starts-with(@name, 'fgd1')]/amplicon/snp">
                            <xsl:if test="snp_call/@count &gt;= $mutant_count_filter and snp_call/@percent &gt;= $prop_filter">
	    		    <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/><xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;]]></xsl:text><xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%),<xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;]]></xsl:text>
                            </xsl:if>
	    		</xsl:for-each></td>
                        <xsl:for-each select="assay[starts-with(@name, 'fgd2')]/amplicon//snp[@name != 'unknown' and @name != 'position of interest']">
		            <td><xsl:choose>
		            <xsl:when test="../significance/@flag"><font color="lightgray"><em><xsl:value-of select="../significance/@flag"/></em></font></xsl:when>
		            <xsl:when test="./significance/@flag"><font color="lightgray"><xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%) - <em><xsl:value-of select="./significance/@flag"/></em></font></xsl:when>
		            <xsl:when test="./significance[not(@flag)]">
		                <xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)
		            </xsl:when>
		            <xsl:otherwise><!-- SNP not present --><em>-</em></xsl:otherwise>
		            </xsl:choose></td>
                        </xsl:for-each>
	    		<td><xsl:for-each select="assay[starts-with(@name, 'fgd2')]/amplicon/snp">
                            <xsl:if test="snp_call/@count &gt;= $mutant_count_filter and snp_call/@percent &gt;= $prop_filter">
	    		    <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/><xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;]]></xsl:text><xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%),<xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;]]></xsl:text>
                            </xsl:if>
	    		</xsl:for-each></td>
                        </tr>
                </xsl:for-each>
            </table><!--</div>-->
	    	<em>Values indicate the number of reads in that sample containing that mutation, out of the total number of reads at that position. A value of '-' indicates that that particular mutation wasn't present in that sample.</em>
	    	<br />
	    	<br />
        </body>
        </html>
    </xsl:template>
    
<!-- Per Sample Detailed Results -->
	<xsl:template match="sample">
        
        <xsl:variable name="prop_filter" select="@proportion_filter * 100"/>
        <xsl:variable name="mutant_count_filter">
            <xsl:choose>
                <xsl:when test="@mutation_depth_filter"><xsl:value-of select="@mutation_depth_filter"/></xsl:when>
                <xsl:otherwise>0</xsl:otherwise>
            </xsl:choose>
        </xsl:variable>

	<exsl:document method="html" href="{/analysis/@run_name}/{@name}.html">
	    <html>
	    <head>
	    	<title>Detailed Results for Sample: <xsl:value-of select="@name"/></title>
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
	        <center><h1>Detailed ASAP Results for Sample: <xsl:value-of select="@name"/></h1></center>
	        <br/>
	    	<h3>SNP profile for sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Assay Name</th>
	    		<th>Average Read Depth</th>
	    		<th>Mutations - # reads containing mutation/depth at position(% reads) - Significance</th>
	    		</tr>
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@type = 'SNP' or @type = 'ROI' or @type = 'mixed'">
	    		    <xsl:call-template name="amplicon-graph"></xsl:call-template>
	    		    <tr>
	    		        <td><a href="#{@name}-graph"><xsl:value-of select="@name"/></a></td>
	    		        <td><xsl:value-of select='format-number(amplicon/average_depth, "#.##")'/></td>
	    		        <xsl:if test="amplicon/@reads &gt; 0">
		    		        <td>
		    		        <xsl:for-each select="amplicon/snp">
		    		            <xsl:if test="./@name != 'unknown'">
                                                <xsl:choose>
                                                  <xsl:when test="./@name = 'position of interest'"><xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/></xsl:when>
                                                  <xsl:otherwise><xsl:value-of select="./@name"/></xsl:otherwise>
                                                </xsl:choose>
                                                - <xsl:value-of select='./snp_call/@count'/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(./snp_call/@percent, "##.##")'/>%)
		    		                <xsl:if test="significance"> - <xsl:value-of select="significance"/><xsl:if test="significance/@flag">(<xsl:value-of select="significance/@flag"/>)</xsl:if></xsl:if>
		    		                <br/>
		    		            </xsl:if>
		    		        </xsl:for-each>
		    		        <xsl:for-each select="amplicon/region_of_interest">
		    		            <xsl:for-each select="mutation">
		    		                <xsl:value-of select="./@name"/> - <xsl:value-of select='./@count'/>(<xsl:value-of select='format-number(./@percent, "##.##")'/>%)
	    		                    <xsl:if test="../significance and ./@percent &gt; $prop_filter"> - <xsl:value-of select="../significance"/><xsl:if test="../significance/@flag">(<xsl:value-of select="../significance/@flag"/>)</xsl:if></xsl:if>
		    		                <br/>
		    		            </xsl:for-each>
		    		        </xsl:for-each>
		    		        </td>
	    		        </xsl:if>
	    		    </tr>
	    		    </xsl:if>
	    		</xsl:for-each>
	    	</table>
	    	<br />
	    	<br />
	    	<a href="../{/analysis/@run_name}.html">Click here to return to the run summary</a>
	    </body>
	    </html>
	</exsl:document>
    </xsl:template>
    
    <xsl:template name="amplicon-graph">
        <div id="{@name}-graph" class="ampGraph">
            <div>
                <a href="#close" title="Close" class="close">X</a>
                <h2>Amplicon Graph</h2>
			<canvas id="{@name}-canvas" height="90vh" class="ampCanvas"></canvas>
			<script>
				var ctx_<xsl:value-of select="translate(translate(@name, '+', '_'), '-', '_')"/> = document.getElementById("<xsl:value-of select="@name"/>-canvas").getContext("2d");
				var chart_<xsl:value-of select="translate(translate(@name, '+', '_'), '-', '_')"/> = new Chart(ctx_<xsl:value-of select="translate(translate(@name, '+', '_'), '-', '_')"/>, {
                                    type: 'bar',
				    data: {
				        labels: "<xsl:value-of select="amplicon/consensus_sequence"/>".split(""),
				        datasets: [{
					    type: 'line',
                                            label: 'Consensus Proportion',
				            yAxisID: 'proportion',
				            data: [<xsl:value-of select="amplicon/proportions"/>],
				            borderColor: "#5F9EA0",
				            borderWidth: 5,
				            fill: false,
				            pointRadius: 0,
				            pointHoverRadius: 3,
				            pointHoverBorderColor: "#B22222",
				        },
				        {
					    type: 'bar',
				            label: 'Read Depth',
				            yAxisID: 'depth',
				            data: [<xsl:value-of select="amplicon/depths"/>],
				            backgroundColor: "#FFDEAD",
				            borderColor: "#DEB887",
				            borderWidth: 1,
				            hoverBorderColor: "#B22222",
				        }]
				    },
				    options: {
					responsive: true,
					hover: {
					    mode: 'label'
					},
					tooltips: {
					    mode: 'label'
					},
				        scales: {
				            yAxes: [{
				            	id: "depth",
				            	position: "left",
				                ticks: {
				                    beginAtZero: true
				                },
				            },
				            {
				            	id: "proportion",
				            	position: "right",
				                ticks: {
				                    beginAtZero: true
				                },
				            }],
				            xAxes: [{
				            	gridLines: {
					            display: false
				            	},
			                        categoryPercentage: 1.0,
			                    }]
				        }
				    }
				});
		    </script>
            </div>
        </div>
    </xsl:template>

</xsl:stylesheet>
