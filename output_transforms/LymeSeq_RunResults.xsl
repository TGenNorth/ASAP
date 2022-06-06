<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:exsl="http://exslt.org/common" xmlns:str="http://exslt.org/strings" xmlns:re="http://exslt.org/regular-expressions" extension-element-prefixes="exsl str re">
    <xsl:import href="http://exslt.org/str/functions/replace/str.replace.function.xsl"/>
    <xsl:output method="html" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" omit-xml-declaration="yes" encoding="UTF-8" indent="yes"/>
    <xsl:template match="/analysis">
        <html>
        <head>
            <title>Run Summary for: <xsl:value-of select="@run_name"/></title>
        </head>
        <body>
            <center><h1>LymeSeq ASAP Run Summary for: <xsl:value-of select="@run_name"/></h1></center>
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
                    <td nowrap="true"><a href="{/analysis/@run_name}/{./@name}.html"><xsl:value-of select="@name"/></a></td>
		    <xsl:for-each select="assay">
		        <td align="center">
			<xsl:choose>
			<xsl:when test="@type = 'gene variant'">
                            <xsl:for-each select="amplicon">
                                <xsl:sort select="@reads" data-type="number" order="descending"/>
                                <xsl:if test="position()=1">
                                    <xsl:value-of select="@variant"/> - <strong><xsl:value-of select="@reads"/></strong>
                                </xsl:if>
                            </xsl:for-each>
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
                            <xsl:variable name="best_amp">
                            <xsl:for-each select="amplicon">
                                <xsl:sort select="breadth" data-type="number" order="descending"/>
                                <xsl:sort select="@reads" data-type="number" order="descending"/>
                                <xsl:if test="position()=1"><xsl:copy-of select="./*"/></xsl:if>
                            </xsl:for-each>
                            </xsl:variable>
			<script>
			function render_<xsl:value-of select="str:replace(str:replace(@name, '+', '_'), '-', '_')"/>() {	
			  var ctx_<xsl:value-of select="str:replace(str:replace(@name, '+', '_'), '-', '_')"/> = document.getElementById("<xsl:value-of select="@name"/>-canvas").getContext("2d");
				var chart_<xsl:value-of select="str:replace(str:replace(@name, '+', '_'), '-', '_')"/> = new Chart(ctx_<xsl:value-of select="str:replace(str:replace(@name, '+', '_'), '-', '_')"/>, {
                                    type: 'bar',
				    data: {
				        labels: "<xsl:value-of select="exsl:node-set($best_amp)/consensus_sequence"/>".split(""),
				        datasets: [{
					    type: 'line',
                                            label: 'Consensus Proportion',
				            yAxisID: 'proportion',
				            data: [<xsl:value-of select="exsl:node-set($best_amp)/proportions"/>],
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
				            data: [<xsl:value-of select="exsl:node-set($best_amp)/depths"/>],
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
	        <center><h1>LymeSeq ASAP Results for Sample: <xsl:value-of select="@name"/></h1></center>
                <table border="0">
                <tr>
                <th>Alignment statistics</th>
                <th>Analysis parameters</th>
                </tr>
                <tr>
	        <td style="padding-left: 20">Total reads: <xsl:value-of select="@mapped_reads + @unmapped_reads + @unassigned_reads"/></td>
                <td style="padding-left: 20">Depth filter: <xsl:value-of select="@depth_filter"/>x</td>
                </tr>
                <tr>
	        <td style="padding-left: 20">Mapped reads: <xsl:value-of select="@mapped_reads"/></td>
                <td style="padding-left: 20">Breadth filter: <xsl:value-of select="@breadth_filter * 100"/>%</td>
                </tr>
                <tr>
	        <td style="padding-left: 20">Unmapped reads: <xsl:value-of select="@unmapped_reads + @unassigned_reads"/></td>
                <td style="padding-left: 20">Proportion filter: <xsl:value-of select="@proportion_filter * 100"/>%</td>
                </tr>
                <tr>
                <td style="padding-left: 20">Aligner used: bowtie2</td>
                </tr>
                </table>
	        <br/>
    	    <xsl:variable name="pre_map">
    	    <xsl:for-each select=".//significance">
                <xsl:if test="not(./@flag)">
                    <entry key="{.}"><xsl:value-of select="ancestor::assay/@name"/></entry>
                </xsl:if>
            </xsl:for-each>
            </xsl:variable>
    	    <xsl:variable name="clinical_map">
    	    <xsl:for-each select="exsl:node-set($pre_map)/entry[not(@key=preceding-sibling::entry/@key)]">
                <xsl:variable name="message" select="@key"/>
                <result message="{$message}">
                <xsl:for-each select="exsl:node-set($pre_map)/entry[@key=$message]">
                    <evidence value="{.}"/>
                </xsl:for-each>
                </result>
            </xsl:for-each>
            </xsl:variable>
	    	<h3>Assays positive for sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" rules="rows" cellpadding="5">
                <xsl:for-each select="exsl:node-set($clinical_map)/result">
    		    <tr>
    		        <td><xsl:text disable-output-escaping="yes"><![CDATA[&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;]]></xsl:text></td>
    		        <td><xsl:value-of select="@message"/> <!--(<xsl:for-each select="evidence[not(@value=preceding-sibling::evidence/@value)]"><xsl:value-of select="@value"/><xsl:if test="position() != last()"><xsl:text>, </xsl:text></xsl:if></xsl:for-each>)--></td>
    		    </tr>
    		    </xsl:for-each>
	    	</table>
	    	<br />
	    	<xsl:variable name="count_16s">
                    <xsl:value-of select="count(./assay[contains(@name, '16S-')]//amplicon[contains(significance, 'burgdorferi')])"/>
                </xsl:variable>
	    	<h4>16S assays positive for Lyme Borreliosis group: <xsl:choose><xsl:when test="$count_16s > 0"><xsl:value-of select="$count_16s+1"/></xsl:when><xsl:otherwise><xsl:value-of select="$count_16s"/></xsl:otherwise></xsl:choose>/4</h4>
	    	<h4>Other assays positive for Lyme Borreliosis group: <xsl:value-of select="count(./assay[not(contains(@name, '16S-'))]/amplicon[contains(significance, 'burgdorferi') or contains(significance, 'Lyme')]/significance[not(@flag)])"/>/5</h4>
	    	<br />
	    	<br />
	    	<h3>Species and strain identification assays for sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Assay Name</th>
	    		<th>Average Read Depth</th>
	    		<th>Coverage Breadth</th>
	    		<th>Significance</th>
	    		<th>SNPs found</th>
	    		</tr>
                        <tbody valign="top">
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="(@function = 'species ID' or @function = 'strain ID') and .//significance/text()"><!--and (exsl:node-set($best_amp)/@reads &gt; 0)--> 
                            <xsl:variable name="best_amp">
                            <xsl:for-each select="amplicon">
                                <xsl:sort select="breadth" data-type="number" order="descending"/>
                                <xsl:sort select="@reads" data-type="number" order="descending"/>
                                <xsl:if test="position()=1"><xsl:copy-of select="./*"/></xsl:if>
                            </xsl:for-each>
                            </xsl:variable>
	    		    <xsl:call-template name="amplicon-graph"></xsl:call-template>
	    		    <tr>
	    		        <td><a href="#{@name}-graph" onclick="render_{str:replace(str:replace(@name, '+', '_'), '-', '_')}()"><xsl:value-of select="@name"/></a></td>
	    		        <td><xsl:value-of select='format-number(exsl:node-set($best_amp)/average_depth, "#.##")'/></td>
	    		        <td><xsl:value-of select='format-number(exsl:node-set($best_amp)/breadth, "##.##")'/>%</td>
	    		        <xsl:choose>
	    		            <xsl:when test="@type = 'presence/absence'">
	    		                <td><xsl:value-of select="amplicon/significance"/><xsl:if test="amplicon/significance/@flag"> (<xsl:value-of select="amplicon/significance/@flag"/>)</xsl:if></td>
	    		            </xsl:when>
	    		            <xsl:when test="@type = 'SNP'">
			    		        <td><xsl:for-each select="amplicon/snp">
			    		            <xsl:if test="significance and ./@name != 'unknown'">
			    		                <xsl:value-of select="./@position"/><xsl:value-of select="./@reference"/>-><xsl:value-of select="./snp_call"/>(<xsl:value-of select='format-number(./snp_call/@percent, "##.##")'/>%)
			    		                - <xsl:value-of select="significance"/><xsl:if test="significance/@flag">(<xsl:value-of select="significance/@flag"/>)</xsl:if>
			    		                <br/>
			    		            </xsl:if>
			    		        </xsl:for-each></td>
	    		            </xsl:when>
	    		            <xsl:when test="@type = 'ROI'">
			    		        <td><xsl:for-each select="amplicon/region_of_interest">
			    		            <xsl:if test="significance">
			    		            <xsl:value-of select="@region"/> - <xsl:value-of select="@reference"/> - <xsl:value-of select="amino_acid_sequence"/>(<xsl:value-of select='format-number(amino_acid_sequence/@percent, "##.##")'/>%)
		    		                - <xsl:value-of select="significance"/>(<xsl:value-of select="significance/@changes"/>)
			    		            <br/>
			    		            </xsl:if>
			    		        </xsl:for-each></td>
	    		            </xsl:when>
                                    <xsl:when test="@type = 'gene variant'">
                                        <td><xsl:for-each select="amplicon">
                                            <xsl:sort select="breadth" data-type="number" order="descending"/>
                                            <xsl:sort select="@reads" data-type="number" order="descending"/>
	    		                    <xsl:if test="breadth &gt; 0">
	    		                        <xsl:value-of select="@variant"/> - <strong><xsl:value-of select="@reads"/></strong> - <xsl:value-of select='format-number(breadth, "##.##")'/>%
	    		                    <xsl:if test="significance">
	    		                        - <xsl:value-of select="significance"/><xsl:if test="significance/@flag"> (<xsl:value-of select="significance/@flag"/>)</xsl:if>
	    		                    </xsl:if>
                                            <br />
		    		            </xsl:if>
		    		        </xsl:for-each></td>
                                    </xsl:when>
                                    <xsl:otherwise><td></td></xsl:otherwise>
	    		        </xsl:choose>
	    		        <td><xsl:if test="exsl:node-set($best_amp)/snp"><div class="snpCell">details...<br/><xsl:for-each select="exsl:node-set($best_amp)/snp">
	    		            <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)<br/>
	    		        </xsl:for-each></div></xsl:if></td>
	    		    </tr>
	    		    </xsl:if>
	    		</xsl:for-each>
                        </tbody>
	    	</table>
	    	<br />
<!--	    	<h3>Antibiotic resistance assays for sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Assay Name</th>
	    		<th>Average Read Depth</th>
	    		<th>Coverage Breadth</th>
	    		<th>Significance</th>
	    		<th>SNPs found</th>
	    		</tr>
                        <tbody valign="top">
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@function = 'resistance type' and .//significance/text()">
                            <xsl:variable name="best_amp">
                            <xsl:for-each select="amplicon">
                                <xsl:sort select="breadth" data-type="number" order="descending"/>
                                <xsl:sort select="@reads" data-type="number" order="descending"/>
                                <xsl:if test="position()=1"><xsl:copy-of select="./*"/></xsl:if>
                            </xsl:for-each>
                            </xsl:variable>
	    		    <xsl:call-template name="amplicon-graph"></xsl:call-template>
	    		    <tr>
	    		        <td><a href="#{@name}-graph" onclick="render_{str:replace(str:replace(@name, '+', '_'), '-', '_')}()"><xsl:value-of select="@name"/></a></td>
	    		        <td><xsl:value-of select='format-number(exsl:node-set($best_amp)/average_depth, "#.##")'/></td>
	    		        <td><xsl:value-of select='format-number(exsl:node-set($best_amp)/breadth, "##.##")'/>%</td>
	    		        <xsl:choose>
	    		            <xsl:when test="@type = 'presence/absence'">
	    		                <td><xsl:value-of select="amplicon/significance"/><xsl:if test="amplicon/significance/@flag"> (<xsl:value-of select="amplicon/significance/@flag"/>)</xsl:if></td>
	    		            </xsl:when>
	    		            <xsl:when test="@type = 'SNP'">
			    		        <td><xsl:for-each select="amplicon/snp">
			    		            <xsl:if test="significance and ./@name != 'unknown'">
			    		                <xsl:value-of select="./@position"/><xsl:value-of select="./@reference"/>-><xsl:value-of select="./snp_call"/>(<xsl:value-of select='format-number(./snp_call/@percent, "##.##")'/>%)
			    		                - <xsl:value-of select="significance"/><xsl:if test="significance/@flag">(<xsl:value-of select="significance/@flag"/>)</xsl:if>
			    		                <br/>
			    		            </xsl:if>
			    		        </xsl:for-each></td>
	    		            </xsl:when>
	    		            <xsl:when test="@type = 'ROI'">
			    		        <td><xsl:for-each select="amplicon/region_of_interest">
			    		            <xsl:if test="significance">
			    		            <xsl:value-of select="@region"/> - <xsl:value-of select="@reference"/> - <xsl:value-of select="amino_acid_sequence"/>(<xsl:value-of select='format-number(amino_acid_sequence/@percent, "##.##")'/>%)
		    		                - <xsl:value-of select="significance"/>(<xsl:value-of select="significance/@changes"/>)
			    		            <br/>
			    		            </xsl:if>
			    		        </xsl:for-each></td>
	    		            </xsl:when>
                                    <xsl:when test="@type = 'gene variant'">
                                        <td><xsl:for-each select="amplicon">
                                            <xsl:sort select="breadth" data-type="number" order="descending"/>
                                            <xsl:sort select="@reads" data-type="number" order="descending"/>
	    		                    <xsl:if test="breadth &gt; 0">
	    		                        <xsl:value-of select="@variant"/> - <strong><xsl:value-of select="@reads"/></strong> - <xsl:value-of select='format-number(breadth, "##.##")'/>%
	    		                    <xsl:if test="significance">
	    		                        - <xsl:value-of select="significance"/><xsl:if test="significance/@flag"> (<xsl:value-of select="significance/@flag"/>)</xsl:if>
	    		                    </xsl:if>
                                            <br />
		    		            </xsl:if>
		    		        </xsl:for-each></td>
                                    </xsl:when>
                                    <xsl:otherwise><td></td></xsl:otherwise>
	    		        </xsl:choose>
	    		        <td><xsl:if test="exsl:node-set($best_amp)/snp"><div class="snpCell">details...<br/><xsl:for-each select="exsl:node-set($best_amp)/snp">
	    		            <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)<br/>
	    		        </xsl:for-each></div></xsl:if></td>
	    		    </tr>
	    		    </xsl:if>
	    		</xsl:for-each>
                        </tbody>
	    	</table>
	    	<br />
	    	<h3>Virulence assays for sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Assay Name</th>
	    		<th>Average Read Depth</th>
	    		<th>Coverage Breadth</th>
	    		<th>Significance</th>
	    		<th>SNPs found</th>
	    		</tr>
                        <tbody valign="top">
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@function = 'virulence factor' and .//significance/text()">
                            <xsl:variable name="best_amp">
                            <xsl:for-each select="amplicon">
                                <xsl:sort select="breadth" data-type="number" order="descending"/>
                                <xsl:sort select="@reads" data-type="number" order="descending"/>
                                <xsl:if test="position()=1"><xsl:copy-of select="./*"/></xsl:if>
                            </xsl:for-each>
                            </xsl:variable>
	    		    <xsl:call-template name="amplicon-graph"></xsl:call-template>
	    		    <tr>
	    		        <td><a href="#{@name}-graph" onclick="render_{str:replace(str:replace(@name, '+', '_'), '-', '_')}()"><xsl:value-of select="@name"/></a></td>
	    		        <td><xsl:value-of select='format-number(exsl:node-set($best_amp)/average_depth, "#.##")'/></td>
	    		        <td><xsl:value-of select='format-number(exsl:node-set($best_amp)/breadth, "##.##")'/>%</td>
	    		        <xsl:choose>
	    		            <xsl:when test="@type = 'presence/absence'">
	    		                <td><xsl:value-of select="amplicon/significance"/><xsl:if test="amplicon/significance/@flag"> (<xsl:value-of select="amplicon/significance/@flag"/>)</xsl:if></td>
	    		            </xsl:when>
	    		            <xsl:when test="@type = 'SNP'">
			    		        <td><xsl:for-each select="amplicon/snp">
			    		            <xsl:if test="significance and ./@name != 'unknown'">
			    		                <xsl:value-of select="./@position"/><xsl:value-of select="./@reference"/>-><xsl:value-of select="./snp_call"/>(<xsl:value-of select='format-number(./snp_call/@percent, "##.##")'/>%)
			    		                - <xsl:value-of select="significance"/><xsl:if test="significance/@flag">(<xsl:value-of select="significance/@flag"/>)</xsl:if>
			    		                <br/>
			    		            </xsl:if>
			    		        </xsl:for-each></td>
	    		            </xsl:when>
	    		            <xsl:when test="@type = 'ROI'">
			    		        <td><xsl:for-each select="amplicon/region_of_interest">
			    		            <xsl:if test="significance">
			    		            <xsl:value-of select="@region"/> - <xsl:value-of select="@reference"/> - <xsl:value-of select="amino_acid_sequence"/>(<xsl:value-of select='format-number(amino_acid_sequence/@percent, "##.##")'/>%)
		    		                - <xsl:value-of select="significance"/>(<xsl:value-of select="significance/@changes"/>)
			    		            <br/>
			    		            </xsl:if>
			    		        </xsl:for-each></td>
	    		            </xsl:when>
                                    <xsl:when test="@type = 'gene variant'">
                                        <td><xsl:for-each select="amplicon">
                                            <xsl:sort select="breadth" data-type="number" order="descending"/>
                                            <xsl:sort select="@reads" data-type="number" order="descending"/>
	    		                    <xsl:if test="breadth &gt; 0">
	    		                        <xsl:value-of select="@variant"/> - <strong><xsl:value-of select="@reads"/></strong> - <xsl:value-of select='format-number(breadth, "##.##")'/>%
	    		                    <xsl:if test="significance">
	    		                        - <xsl:value-of select="significance"/><xsl:if test="significance/@flag"> (<xsl:value-of select="significance/@flag"/>)</xsl:if>
	    		                    </xsl:if>
                                            <br />
		    		            </xsl:if>
		    		        </xsl:for-each></td>
                                    </xsl:when>
                                    <xsl:otherwise><td></td></xsl:otherwise>
	    		        </xsl:choose>
	    		        <td><xsl:if test="exsl:node-set($best_amp)/snp"><div class="snpCell">details...<br/><xsl:for-each select="exsl:node-set($best_amp)/snp">
	    		            <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)<br/>
	    		        </xsl:for-each></div></xsl:if></td>
	    		    </tr>
	    		    </xsl:if>
	    		</xsl:for-each>
                        </tbody>
	    	</table>
-->                <br />
                <center><img src="../web_resources/TGen-North_small_logo.png"/></center>
	    </body>
	    </html>
	</exsl:document>
    </xsl:template>
</xsl:stylesheet>
