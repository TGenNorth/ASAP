<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:asap="http://pathogen.tgen.org/ASAP/functions" xmlns:exsl="http://exslt.org/common" xmlns:str="http://exslt.org/strings" extension-element-prefixes="exsl str asap">
    <xsl:output method="xhtml" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" omit-xml-declaration="yes" encoding="UTF-8" indent="yes"/>

<!-- Clinical Run Summary -->
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
        <xsl:variable name="prop_filter" select="sample[1]/@proportion_filter * 100"/>
        <xsl:variable name="mutant_count_filter">
            <xsl:choose>
                <xsl:when test="sample[1]/@mutation_depth_filter"><xsl:value-of select="sample[1]/@mutation_depth_filter"/></xsl:when>
                <xsl:otherwise>0</xsl:otherwise>
            </xsl:choose>
        </xsl:variable>
        	<center><h1>ASAP Run Summary for: <xsl:value-of select="@run_name"/></h1></center>
	        <br />
	        <br />
            <div class="div-table-column-locked"><table class="table-column-locked">
	    		<tr>
	    		<th class="headcol">Mutation</th>
	    		<xsl:for-each select="sample">
	    		    <th nowrap="true">
	    		        <a href="{/analysis/@run_name}/{./@name}.html"><xsl:value-of select="@name"/></a>
	    		    </th>
	    		</xsl:for-each>
	    		</tr>
	    		<xsl:for-each select="sample[1]/assay[not(@type='presence/absence')]">
	    		    <xsl:variable name="current_assay" select="./@name"/>
                            <xsl:text disable-output-escaping="yes"><![CDATA[<tr>]]></xsl:text>
	    		    <xsl:if test="@type = 'SNP' or @type = 'mixed'">
                            <xsl:for-each select="asap:distinct-values(//assay[@name=current()/@name]//amplicon//snp/@name[. != 'unknown' and . != 'position of interest'])">
                            <xsl:sort select="."/>
	    		            <xsl:variable name="current_snp" select="."/>
	                            <th class="headcol"><xsl:value-of select="."/></th>
	                            <xsl:for-each select="//sample">
	                                <td nowrap="true">
	                                <xsl:if test="not(.//assay[@name=$current_assay]//amplicon//snp[@name=$current_snp])"><!-- assay not present --><em>no coverage</em></xsl:if>
	                                <xsl:for-each select=".//assay[@name=$current_assay]//amplicon//snp[@name=$current_snp]">
                                            <xsl:if test="position()=1">
		                            <xsl:choose>
		                            <xsl:when test="../significance/@flag"><em><xsl:value-of select="../significance/@flag"/></em></xsl:when>
		                            <xsl:when test="./significance/@flag"><em><xsl:value-of select="./significance/@flag"/></em></xsl:when>
		                            <xsl:otherwise>
                                                <xsl:variable name="DEPTH" select="@depth"/>
                                                <xsl:variable name="REF" select="@reference"/>
                                                <xsl:for-each select="base_distribution/@*">
                                                <xsl:if test="name() = $REF"><xsl:text disable-output-escaping="yes"><![CDATA[<em>]]></xsl:text></xsl:if>
                                                    <xsl:value-of select="concat(name(), ':', .,'/')"/><xsl:value-of select="$DEPTH"/>(<xsl:value-of select='format-number(. div $DEPTH * 100, "##.###")'/>%)
                                                <xsl:if test="name() = $REF"><xsl:text disable-output-escaping="yes"><![CDATA[</em>]]></xsl:text></xsl:if>
                                                </xsl:for-each>
                                            </xsl:otherwise>
		                            </xsl:choose>
                                            </xsl:if>
	                                </xsl:for-each>
	                                </td>
	                            </xsl:for-each>
	                            <xsl:text disable-output-escaping="yes"><![CDATA[</tr><tr>]]></xsl:text>
                            </xsl:for-each>
	    		    </xsl:if>
	    		    <xsl:if test="@type = 'ROI' or @type = 'mixed'">
                            <xsl:for-each select="asap:distinct-values(//assay[@name=current()/@name]//amplicon//region_of_interest//mutation/@name)">
                            <xsl:sort select="."/>
	    		            <xsl:variable name="current_codon" select="."/>
	                            <th class="headcol"><xsl:value-of select="$current_assay"/>-SMOR <xsl:value-of select="."/></th>
	                            <xsl:for-each select="//sample">
	                                <td nowrap="true">
	                                <xsl:if test="not(.//assay[@name=$current_assay]//amplicon//region_of_interest//mutation[@name=$current_codon])"><!-- assay not present --><em>no coverage</em></xsl:if>
	                                <xsl:for-each select=".//assay[@name=$current_assay]//amplicon//region_of_interest//mutation[@name=$current_codon]">
		                            <xsl:choose>
		                            <xsl:when test="../../significance/@flag"><em><xsl:value-of select="../../significance/@flag"/></em></xsl:when>
		                            <xsl:when test="../significance/@flag"><em><xsl:value-of select="../significance/@flag"/></em></xsl:when>
		                            <xsl:when test="../significance[not(@flag)] and @count &gt;= $mutant_count_filter and @percent &gt;= $prop_filter">
		                                <xsl:value-of select="@count"/>/<xsl:value-of select="../@depth"/>(<xsl:value-of select='format-number(@percent, "##.##")'/>%)
		                            </xsl:when>
		                            <xsl:otherwise><!-- mutant codon not present --><xsl:value-of select="@count"/>/<xsl:value-of select="../@depth"/>(<xsl:value-of select='format-number(@percent, "##.##")'/>%)</xsl:otherwise>
		                            </xsl:choose>
	                                </xsl:for-each>
	                                </td>
	                            </xsl:for-each>
                            <xsl:text disable-output-escaping="yes"><![CDATA[</tr><tr>]]></xsl:text>
                            </xsl:for-each>
	    		    </xsl:if>
	    		    <xsl:text disable-output-escaping="yes"><![CDATA[</tr>]]></xsl:text>
	    		<!--<xsl:apply-templates select="."/>  -->
                </xsl:for-each>
            </table></div>
	    	<em>Values indicate the number of reads in that sample containing that mutation.</em>
	    	<br />
	    	<br />
	    	<a href="{@run_name}_details.html">Click here for more details</a>
        </body>
        </html>
    
<!-- Detailed Run Summary -->
    <exsl:document method="html" href="{@run_name}_details.html">
        <html>
        <head>
            <title>Run Summary for: <xsl:value-of select="@run_name"/></title>
        </head>
        <body>
        	<center><h1>ASAP Run Summary for: <xsl:value-of select="@run_name"/></h1></center>
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
	    	<br />
	    	<br />
	    	<a href="{@run_name}.html">Click here for clinical summary</a>
        </body>
        </html>
    </exsl:document>
	</xsl:template>
	
	<xsl:template match="sample">
        
        <xsl:variable name="prop_filter" select="@proportion_filter * 100"/>
        <xsl:variable name="mutant_count_filter">
            <xsl:choose>
                <xsl:when test="@mutation_depth_filter"><xsl:value-of select="@mutation_depth_filter"/></xsl:when>
                <xsl:otherwise>0</xsl:otherwise>
            </xsl:choose>
        </xsl:variable>

<!-- Per Sample Clinical Results -->
	<exsl:document method="html" href="{/analysis/@run_name}/{@name}.html">
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
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Rifampin') and contains(@level, 'high')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Rifampin') and not(@level='low')]"><font color="red">HR</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Rifampin')]"><font color="red">LHR</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Rifampin')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Isoniazid') and contains(@level, 'high')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Isoniazid') and not(@level='low')]"><font color="red">HR</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Isoniazid')]"><font color="red">LHR</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Isoniazid')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Quinolones') and contains(@level, 'high')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Quinolones') and not(@level='low')]"><font color="red">HR</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Quinolones')]"><font color="red">LHR</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Quinolones')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Kanamycin') and contains(@level, 'high')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Kanamycin') and not(@level='low')]"><font color="red">HR</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Kanamycin')]"><font color="red">LHR</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Kanamycin')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Capreomycin') and contains(@level, 'high')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Capreomycin') and not(@level='low')]"><font color="red">HR</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Capreomycin')]"><font color="red">LHR</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Capreomycin')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Amikacin') and contains(@level, 'high')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Amikacin') and not(@level='low')]"><font color="red">HR</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Amikacin')]"><font color="red">LHR</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Amikacin')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Pyrazinamide') and contains(@level, 'high')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Pyrazinamide') and not(@level='low')]"><font color="red">HR</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Pyrazinamide')]"><font color="red">LHR</font></xsl:when>
                    	<xsl:when test=".//significance[@flag and contains(@resistance, 'Pyrazinamide')]">Ind.</xsl:when>
                    	<xsl:otherwise>S</xsl:otherwise>
                    </xsl:choose></td>
                    <td class="rotate" align="center"><xsl:choose>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Ethambutal') and contains(@level, 'high')]"><font color="red">R</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Ethambutal') and not(@level='low')]"><font color="red">HR</font></xsl:when>
                    	<xsl:when test=".//significance[not(@flag) and contains(@resistance, 'Ethambutal')]"><font color="red">LHR</font></xsl:when>
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
	    		        <xsl:when test="@type = 'mixed'"><th><xsl:value-of select="@name"/> SNP/codon (% res)</th></xsl:when>
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
	    		                    <xsl:if test="significance">
	    		                        <xsl:for-each select="mutation">
	    		                            <xsl:if test="@percent &gt; $prop_filter and @count &gt;= $mutant_count_filter"><xsl:value-of select="@name"/> (<xsl:value-of select='format-number(@percent, "##.##")'/>%)<br/></xsl:if>
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
	</exsl:document>

<!-- Per Sample Detailed Results -->
	<exsl:document method="html" href="{/analysis/@run_name}/{@name}_details.html">
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
	    		    <tr>
	    		        <td><xsl:value-of select="@name"/></td>
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
	    		    <tr>
	    		        <td><xsl:value-of select="@name"/></td>
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
    
</xsl:stylesheet>
