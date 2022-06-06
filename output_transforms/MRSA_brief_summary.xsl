<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" xmlns:str="http://exslt.org/strings" extension-element-prefixes="exsl str">
    <xsl:import href="http://exslt.org/str/functions/replace/str.replace.function.xsl"/>
    <xsl:output method="html" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" omit-xml-declaration="yes" encoding="UTF-8" indent="yes"/>
    <xsl:template match="/analysis">    
        <html>
        <head>
            <title>Run Summary for: <xsl:value-of select="@run_name"/></title>
            <style type="text/css">
				.table-header-rotated {
				  border-collapse: collapse;
				}
				.table-header-rotated td.rotate {
				  width: 45px;
				}
				.table-header-rotated td.norotate {
				  text-align: left;
				  white-space: nowrap;
				}
				.table-header-rotated th.norotate {
				  padding: 10px 40px;
				  vertical-align: bottom;
				}
				.table-header-rotated td {
				  text-align: center;
				  padding: 10px 5px;
				  border: 2px solid #aaa;
				}
				.table-header-rotated th.rotate {
				  height: 140px;
				  white-space: nowrap;
				}
				.table-header-rotated th.rotate <xsl:text disable-output-escaping="yes"><![CDATA[>]]></xsl:text> div {
				  -webkit-transform: translate(40px, 51px) rotate(315deg);
				      -ms-transform: translate(40px, 51px) rotate(315deg);
				          transform: translate(40px, 51px) rotate(315deg);
				  width: 30px;
				}
				.table-header-rotated th.rotate <xsl:text disable-output-escaping="yes"><![CDATA[>]]></xsl:text> div <xsl:text disable-output-escaping="yes"><![CDATA[>]]></xsl:text> span {
				  border-bottom: 2px solid #aaa;
				  padding: 5px 10px;
				}
				.table-header-rotated th.row-header {
				  padding: 0 10px;
				  border-bottom: 2px solid #aaa;
                                }
            </style>
        </head>
        <body>
        	<center><h1>ASAP Run Summary for: <xsl:value-of select="@run_name"/></h1></center>
	        <br/>
            <table class="table-header-rotated">
	    		<tr>
	    		<th class="norotate">Sample</th>
	    		<th class="rotate"><div><span><em>Staphylococcus aureus</em> Confirmed</span></div></th>
	    		<th class="rotate"><div><span>MRSA Confirmed</span></div></th>
	    		<th class="norotate">Lineage</th>
	    		<th class="rotate"><div><span>Aminoglycoside Resistance</span></div></th>
	    		<th class="rotate"><div><span>β-lactam Resistance</span></div></th>
	    		<th class="rotate"><div><span>Clindamycin Resistance</span></div></th>
	    		<th class="rotate"><div><span>Fluoroquinolone Resistance</span></div></th>
	    		<th class="rotate"><div><span>Linezolid Resistance</span></div></th>
	    		<th class="rotate"><div><span>Macrolide Resistance</span></div></th>
	    		<th class="rotate"><div><span>Rifampin Resistance</span></div></th>
	    		<th class="rotate"><div><span>Tetracycline Resistance</span></div></th>
	    		<th class="rotate"><div><span>Vancomycin Resistance</span></div></th>
	    		</tr>
                <xsl:for-each select="sample">
                    <tr>
                        <td class="norotate"><a href="{/analysis/@run_name}/{./@name}.html"><xsl:value-of select="@name"/></a></td>
                        <td class="rotate" align="center"><xsl:choose><xsl:when test=".//significance[not(@flag)]='Staphylococcus aureus present'"><img src="web_resources/check.png" style="width:30px;height:30px;"/></xsl:when><xsl:otherwise><img src="web_resources/cross.png" style="width:30px;height:30px;"/></xsl:otherwise></xsl:choose></td>
                        <td class="rotate" align="center"><xsl:choose><xsl:when test=".//significance[not(@flag)]='MRSA present'"><img src="web_resources/check.png" style="width:30px;height:30px;"/></xsl:when><xsl:otherwise><img src="web_resources/cross.png" style="width:30px;height:30px;"/></xsl:otherwise></xsl:choose></td>
                        <td class="norotate"><xsl:choose>
                        <xsl:when test=".//significance[not(@flag)]='S. aureus CC5 subtype present'">CC5</xsl:when>
                        <xsl:when test=".//significance[not(@flag) and contains(., 'S. aureus CC8')]">CC8
                            <xsl:choose>
                                <xsl:when test=".//significance[not(@flag)]='S. aureus ST239 subtype present'"> -> ST239</xsl:when>
                                <xsl:when test=".//significance[not(@flag)]='S. aureus Archaic Iberian subtype present'"> -> Archaic Iberian</xsl:when>
                                <xsl:when test=".//significance[not(@flag)]='S. aureus True Iberian subtype present'"> -> True Iberian</xsl:when>
                                <xsl:when test=".//significance[not(@flag)]='S. aureus USA300/500 subtype present'"> -> USA300 or USA500
                                    <xsl:if test=".//significance[not(@flag)]='S. aureus USA300 subtype present'"> -> USA300</xsl:if>                                  
                                </xsl:when>
                            </xsl:choose>
                        </xsl:when>
                        <xsl:otherwise><em>Undetermined</em></xsl:otherwise></xsl:choose></td>
                        <td class="rotate" align="center"><xsl:choose>
                            <xsl:when test=".//significance[not(@flag)]='Aminoglycoside resistance'"><font color="red">R</font></xsl:when>
                            <xsl:when test=".//significance[not(@flag)]='Resistant to all aminoglycosides (e.g. Gentamicin, Amikacin, Tobramycin)'"><font color="red">R</font></xsl:when>
                            <xsl:otherwise>S</xsl:otherwise></xsl:choose></td>
                        <td class="rotate" align="center"><xsl:choose><xsl:when test=".//significance[not(@flag) and contains(., 'β-lactam resistance')]"><font color="red">R</font></xsl:when><xsl:otherwise>S</xsl:otherwise></xsl:choose></td>
                        <td class="rotate" align="center"><xsl:choose><xsl:when test=".//significance[not(@flag)]='Clindamycin resistance'"><font color="red">R</font></xsl:when><xsl:otherwise>S</xsl:otherwise></xsl:choose></td>
                        <td class="rotate" align="center"><xsl:choose><xsl:when test=".//significance[not(@flag)]='Fluoroquinolone resistance'"><font color="red">R</font></xsl:when><xsl:otherwise>S</xsl:otherwise></xsl:choose></td>
                        <td class="rotate" align="center"><xsl:choose><xsl:when test=".//significance[not(@flag)]='Linezolid resistance'"><font color="red">R</font></xsl:when><xsl:otherwise>S</xsl:otherwise></xsl:choose></td>
                        <td class="rotate" align="center"><xsl:choose><xsl:when test=".//significance[not(@flag)]='Macrolide resistance'"><font color="red">R</font></xsl:when><xsl:otherwise>S</xsl:otherwise></xsl:choose></td>
                        <td class="rotate" align="center"><xsl:choose><xsl:when test=".//significance[not(@flag)]='Rifampin resistance'"><font color="red">R</font></xsl:when><xsl:otherwise>S</xsl:otherwise></xsl:choose></td>
                        <td class="rotate" align="center"><xsl:choose><xsl:when test=".//significance[not(@flag)]='Tetracycline resistance'"><font color="red">R</font></xsl:when><xsl:otherwise>S</xsl:otherwise></xsl:choose></td>
                        <td class="rotate" align="center"><xsl:choose>
                            <xsl:when test=".//significance[not(@flag)]='Vancomycin resistance'"><font color="red">R</font></xsl:when>
                            <xsl:when test=".//significance[not(@flag)]='Vancomycin intermediate resistance'"><font color="red">I</font></xsl:when>
                            <xsl:otherwise>S</xsl:otherwise></xsl:choose></td>
                    </tr>
                    <xsl:apply-templates select="."/>
                </xsl:for-each>
            </table>
            <br />
            <center><img src="web_resources/TGen-North_small_logo.png"/></center>
        </body>
        </html>
    </xsl:template>

    <xsl:template name="amplicon-graph">
        <div id="{@name}-graph" class="ampGraph">
            <div>
                <a href="#close" title="Close" class="close">X</a>
                <h2>Amplicon Graph</h2>
			<canvas id="{@name}-canvas" height="90vh" class="ampCanvas"></canvas>
			<script>
			  function render_<xsl:value-of select="str:replace(str:replace(@name, '+', '_'), '-', '_')"/>() {
			  var ctx_<xsl:value-of select="str:replace(str:replace(@name, '+', '_'), '-', '_')"/> = document.getElementById("<xsl:value-of select="@name"/>-canvas").getContext("2d");
				var chart_<xsl:value-of select="str:replace(str:replace(@name, '+', '_'), '-', '_')"/> = new Chart(ctx_<xsl:value-of select="str:replace(str:replace(@name, '+', '_'), '-', '_')"/>, {
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
				}
		    </script>
            </div>
        </div>
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
	    	<h3>Presence/absence assays present in sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Assay Name</th>
	    		<th># of Reads</th>
	    		<th>Coverage Breadth</th>
	    		<th>Significance</th>
	    		<th>SNPs found(% reads containing SNP)</th>
	    		</tr>
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@type = 'presence/absence' and amplicon/@reads &gt; 0">
	    		    <xsl:call-template name="amplicon-graph"></xsl:call-template>
	    		    <tr>
	    		        <td><a href="#{@name}-graph" onclick="render_{str:replace(str:replace(@name, '+', '_'), '-', '_')}()"><xsl:value-of select="@name"/></a></td>
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
	    		    <xsl:call-template name="amplicon-graph"></xsl:call-template>
	    		    <tr>
	    		        <td><a href="#{@name}-graph" onclick="render_{str:replace(str:replace(@name, '+', '_'), '-', '_')}()"><xsl:value-of select="@name"/></a></td>
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
	    	<h3>Region of Interest assays for sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Assay Name</th>
	    		<th># of Reads</th>
	    		<th>Coverage Breadth</th>
	    		<th>Region - Reference - Most common sequence (% reads containing seq) - Significance (# of aa changes)</th>
	    		<th>SNPs found(% reads containing SNP)</th>
	    		</tr>
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@type = 'ROI' or @type = 'mixed' and amplicon/region_of_interest">
	    		    <xsl:call-template name="amplicon-graph"></xsl:call-template>
	    		    <tr>
	    		        <td><a href="#{@name}-graph" onclick="render_{str:replace(str:replace(@name, '+', '_'), '-', '_')}()"><xsl:value-of select="@name"/></a></td>
	    		        <td><xsl:value-of select="amplicon/@reads"/></td>
	    		        <xsl:if test="amplicon/@reads &gt; 0">
		    		        <td><xsl:value-of select='format-number(amplicon/breadth, "##.##")'/>%</td>
		    		        <td><xsl:for-each select="amplicon/region_of_interest">
		    		            <xsl:value-of select="@region"/> - <xsl:value-of select="@reference"/> - <xsl:value-of select="amino_acid_sequence"/>(<xsl:value-of select='format-number(amino_acid_sequence/@percent, "##.##")'/>%)
	    		                - <xsl:value-of select="significance"/>(<xsl:value-of select="significance/@changes"/>)
		    		            <br/>
		    		        </xsl:for-each></td>
		    		        <td><xsl:for-each select="amplicon/snp">
		    		            <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)<br/>
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
	    </body>
	    </html>
	</exsl:document>
    </xsl:template>
</xsl:stylesheet>
