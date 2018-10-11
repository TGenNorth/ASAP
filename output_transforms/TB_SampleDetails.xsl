<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" xmlns:str="http://exslt.org/strings" extension-element-prefixes="exsl str">
    <xsl:output method="xhtml" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" omit-xml-declaration="yes" encoding="UTF-8" indent="yes"/>

	<xsl:template match="sample">

        <xsl:variable name="prop_filter" select="@proportion_filter * 100"/>

<!-- Per Sample Clinical Results -->
	    <html>
	    <head>
	    	<title>Clinical Results for Sample: <xsl:value-of select="@name"/></title>
            <style>
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
	        <center><h1>TB Clinical ASAP Results for Sample: <xsl:value-of select="@name"/></h1></center>
	        <br />
	        <br />
            <table class="table-header-rotated">
	    		<tr>
	    		<th class="norotate">Sample</th>
	    		<th class="rotate"><div><span><em>M. tubercolosis</em> Confirmed</span></div></th>
	    		<th class="rotate"><div><span>Rifampin Resistance</span></div></th>
	    		<th class="rotate"><div><span>Isoniazid Resistance</span></div></th>
	    		<th class="rotate"><div><span>Quinolone Resistance</span></div></th>
	    		<th class="rotate"><div><span>Kanamycin Resistance</span></div></th>
	    		<th class="rotate"><div><span>Capreomycin Resistance</span></div></th>
	    		<th class="rotate"><div><span>Amikacin Resistance</span></div></th>
	    		<th class="rotate"><div><span>Pyrazinamide Resistance</span></div></th>
	    		<th class="rotate"><div><span>Ethambutal Resistance</span></div></th>
	    		</tr>
                <tr>
                    <td class="norotate"><xsl:value-of select="@name"/></td>
                    <td class="rotate" align="center"><xsl:choose><xsl:when test=".//significance[not(@flag)]='Mycobacterium tuberculosis confirmed'"><img src="../check.png" style="width:30px;height:30px;"/></xsl:when><xsl:otherwise><img src="../cross.png" style="width:30px;height:30px;"/></xsl:otherwise></xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Rifampin')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Rifampin')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Isoniazid')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Isoniazid')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Quinolones')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Quinolones')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Kanamycin')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Kanamycin')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Capreomycin')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Capreomycin')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Amikacin')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Amikacin')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Pyrazinamide')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Pyrazinamide')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Ethambutal')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Ethambutal')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                </tr>
            </table>
            <br />
            <br />
	    	<table border="2" rules="rows">
	    	    <tr><th colspan="2">Clinical Summary for Sample: <xsl:value-of select="@name"/></th></tr>
	    	    <xsl:for-each select=".//significance">
                        <xsl:if test="not(./@flag) and not(./@changes)">
	    		    <tr>
	    		    <td><xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;]]></xsl:text></td>
	    		    <td><xsl:value-of select="."/></td>
	    		    </tr>
                        </xsl:if>
	            </xsl:for-each>
	    	</table>
	    	<br />
	    	<br />
	    	<h3>Resistance Mutations Present in Sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Gene Target:</th>
	    		<xsl:for-each select="assay">
	    		    <xsl:choose>
	    		        <xsl:when test="@type = 'SNP'"><th><xsl:value-of select="@name"/> SNP (% res)</th></xsl:when>
	    		        <xsl:when test="@type = 'ROI'"><th><xsl:value-of select="@name"/> codon (% res)</th></xsl:when>
	    		        <xsl:otherwise></xsl:otherwise>
	    		    </xsl:choose>
	    		</xsl:for-each>
	    		</tr>
	    		<tr>
	    		<th>Mutations:</th>
	    		<xsl:for-each select="assay[not(@type='presence/absence')]">
	    		    <td align="center"><xsl:for-each select="amplicon">
	    		        <xsl:choose>
			    	    <xsl:when test=".//significance/@flag"><em><xsl:value-of select=".//significance/@flag"/></em></xsl:when>
	    		            <xsl:when test="snp/significance">
	    		                <xsl:for-each select="snp">
			    		    <xsl:if test="significance">
                                                <xsl:choose>
                                                    <xsl:when test="@name = 'position of interest'"><xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/></xsl:when>
                                                    <xsl:otherwise><xsl:value-of select="@name"/></xsl:otherwise>
                                                </xsl:choose>
                                                (<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)<br/>
                                            </xsl:if>
	    		                </xsl:for-each>
	    		            </xsl:when>
	    		            <xsl:when test="region_of_interest/significance">
	    		                <xsl:for-each select="region_of_interest">
	    		                    <xsl:if test="significance and not(significance/@changes)">
	    		                        <xsl:for-each select="mutation">
	    		                            <xsl:if test="@percent &gt; $prop_filter"><xsl:value-of select="@name"/> (<xsl:value-of select='format-number(@percent, "##.##")'/>%)<br/></xsl:if>
	    		                        </xsl:for-each>
	    		                    </xsl:if>
	    		                </xsl:for-each>
	    		            </xsl:when>
	    		            <xsl:otherwise><em>none</em></xsl:otherwise>
	    		        </xsl:choose>
	    		    </xsl:for-each></td>
	    		</xsl:for-each>
	    		</tr>
	    	</table>
	    	<em>Percentages indicate the percentage of the sample containing that mutation, a value of 'none' indicates that no resistant mutations were present in that gene.</em>
	    	<br />
	    	<br />
	    	<a href="{@name}_details.html">Click here for more sample details</a>
	    	<br />
	    	<br />
	    	<a href="../{/analysis/@run_name}.html">Click here to return to the run summary</a>
	    </body>
	    </html>

<!-- Per Sample Detailed Results -->
	<exsl:document method="html" href="{@name}_details.html">
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
	    		<th>Average Read Depth</th>
	    		<th>Coverage Breadth</th>
	    		<th>Significance</th>
	    		<th>SNPs found(% reads containing SNP)</th>
	    		</tr>
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@name = 'IS6110'">
	    		    <xsl:call-template name="amplicon-graph"></xsl:call-template>
	    		    <tr>
	    		        <td><a href="#{@name}-graph"><xsl:value-of select="@name"/></a></td>
	    		        <td><xsl:value-of select='format-number(amplicon/average_depth, "#.##")'/></td>
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
	    		<th>Average Read Depth</th>
	    		<th>Known Mutations - # reads containing mutation(% reads) - Significance</th>
	    		<th>Other SNPs found(% reads containing SNP)</th>
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
                                                - <xsl:value-of select='./snp_call/@count'/>(<xsl:value-of select='format-number(./snp_call/@percent, "##.##")'/>%)
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
	    	<br />
	    	<a href="{@name}.html">Click here for clinical sample results</a>
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
		    </script>
            </div>
        </div>
    </xsl:template>

</xsl:stylesheet>
