<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:exsl="http://exslt.org/common" xmlns:str="http://exslt.org/strings" extension-element-prefixes="exsl str">
  <xsl:import href="http://exslt.org/str/functions/replace/str.replace.function.xsl"/>
    <xsl:output method="xhtml" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" omit-xml-declaration="yes" encoding="UTF-8" indent="yes"/>
    <xsl:template match="/analysis">
        <html>
        <head>
            <title>Run Summary for: <xsl:value-of select="@run_name"/></title>
        </head>
        <body>
            <center><h1>ASAP Run Summary for: <xsl:value-of select="@run_name"/></h1></center>
	    <br/>
	    <em>Total number of reads aligning</em>
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
                    <td nowrap="true"><a href="{/analysis/@run_name}/{./@name}.html"><xsl:value-of select="@name"/></a></td>
		    <xsl:for-each select="assay">
		        <td align="center">
			<xsl:choose>
			<xsl:when test="@type = 'gene variant'">
                            <xsl:value-of select='sum(./amplicon/@reads)'/>
			</xsl:when>
			<xsl:otherwise><xsl:value-of select="amplicon/@reads"/></xsl:otherwise>
			</xsl:choose>
			</td>
		    </xsl:for-each>
                    </tr>
                    <xsl:apply-templates select="."/>
                </xsl:for-each>
            </table>
        </body>
        </html>
    </xsl:template>

    <xsl:template name="amplicon-graph">
        <div id="{@name}-graph" class="ampGraph">
            <div>
                <a href="#close" title="Close" class="close">X</a>
                <h2><xsl:value-of select="@name"/> Coverage Graph</h2>
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
				            	id: "proportion",
				            	position: "left",
				                ticks: {
				                    beginAtZero: true
				                },
				            },
				            {
				            	id: "depth",
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
		    .snpCell{
                height:15px;
                overflow:hidden;
                text-overflow:ellipsis
            }
            .snpCell:hover{
                height:auto;
                width:auto;
            }
	    	</style>
	    </head>
	    <body>
	        <center><h1>ASAP Results for Sample: <xsl:value-of select="@name"/></h1></center>
                <table border="0">
                <tr>
                <th>Alignment statistics</th>
                <th>Analysis parameters</th>
                </tr>
                <tr>
	        <td style="padding-left: 20">Total reads: <xsl:value-of select="@mapped_reads + @unmapped_reads"/></td>
                <td style="padding-left: 20">Depth filter: <xsl:value-of select="@depth_filter"/>x</td>
                </tr>
                <tr>
	        <td style="padding-left: 20">Mapped reads: <xsl:value-of select="@mapped_reads"/></td>
                <td style="padding-left: 20">Breadth filter: <xsl:value-of select="@breadth_filter * 100"/>%</td>
                </tr>
                <tr>
	        <td style="padding-left: 20">Unmapped reads: <xsl:value-of select="@unmapped_reads"/></td>
                <td style="padding-left: 20">Proportion filter: <xsl:value-of select="@proportion_filter * 100"/>%</td>
                </tr>
                <tr>
                <td style="padding-left: 20">Aligner used: bowtie2</td>
                <td style="padding-left: 20">Identity filter: <xsl:value-of select="@identity_filter * 100"/>%</td>
                </tr>
                </table>
	        <br/>
	    	<table border="2" rules="rows">
	    	    <tr><th colspan="2">Summary for Sample: <xsl:value-of select="@name"/></th></tr>
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
	    	<h3>Details for sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Assay Name</th>
	    		<th># of Reads(# failed identity filter)</th>
	    		<th>Coverage Breadth</th>
	    		<th>Significance</th>
	    		<th>SNPs found(% reads containing SNP)</th>
	    		</tr>
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@type = 'presence/absence'">
	    		    <xsl:call-template name="amplicon-graph"></xsl:call-template>
	    		    <tr>
	    		        <td><a href="#{@name}-graph" onclick="render_{str:replace(str:replace(@name, '+', '_'), '-', '_')}()"><xsl:value-of select="@name"/></a></td>
	    		        <td><xsl:value-of select="amplicon/@reads"/>(<xsl:value-of select="amplicon/@discarded_reads"/>)</td>
	    		        <xsl:choose><xsl:when test="amplicon/breadth">
                                  <td><xsl:value-of select='format-number(amplicon/breadth, "##.##")'/>%</td>
                                </xsl:when><xsl:otherwise>
                                  <td>0%</td>
                                </xsl:otherwise></xsl:choose>
	    		        <td><xsl:value-of select="amplicon/significance"/><xsl:if test="amplicon/significance/@flag"> (<xsl:value-of select="amplicon/significance/@flag"/>)</xsl:if></td>
	    		        <td><xsl:for-each select="amplicon/snp">
	    		            <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)<br/>
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
