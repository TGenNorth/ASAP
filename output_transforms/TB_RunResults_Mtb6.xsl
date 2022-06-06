<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common" extension-element-prefixes="exsl">
    <xsl:output method="xhtml" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" omit-xml-declaration="yes" encoding="UTF-8" indent="yes"/>

<!-- Clinical Run Summary -->
    <xsl:template match="/analysis">
        <html>
        <head>
            <title>Run Summary for: <xsl:value-of select="@run_name"/></title>
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
        	<center><h1>TB Clinical ASAP Run Summary for: <xsl:value-of select="@run_name"/></h1></center>
	        <br />
	        <br />
            <table class="table-header-rotated">
	    		<tr>
	    		<th class="norotate">Sample</th>
	    		<th class="rotate"><div><span>Rifampin Resistance</span></div></th>
	    		<th class="rotate"><div><span>Isoniazid Resistance</span></div></th>
	    		<th class="rotate"><div><span>Quinolone Resistance</span></div></th>
	    		<th class="rotate"><div><span>Kanamycin Resistance</span></div></th>
	    		<th class="rotate"><div><span>Capreomycin Resistance</span></div></th>
	    		<th class="rotate"><div><span>Amikacin Resistance</span></div></th>
	    		</tr>
                <xsl:for-each select="sample">
                    <tr>
                        <td class="norotate"><a href="{/analysis/@run_name}/{./@name}.html"><xsl:value-of select="@name"/></a></td>
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
                    </tr>
                </xsl:for-each>
            </table>
            <br />
            <br />
            <table border="1" cellpadding="3">
	    		<tr>
	    		<th>Sample</th>
	    		<xsl:for-each select="sample[position()&lt;=1]">
	    		<xsl:for-each select="assay[not(@type='presence/absence')]">
	    		    <th nowrap="true">
	    		    <xsl:choose>
	    		        <xsl:when test="@type = 'SNP'"><xsl:value-of select="@name"/> SNP (% res)</xsl:when>
	    		        <xsl:when test="@type = 'ROI'"><xsl:value-of select="@name"/> codon (% res)</xsl:when>
	    		        <xsl:otherwise><xsl:value-of select="@name"/></xsl:otherwise>
	    		    </xsl:choose>
	    		    </th>
	    		</xsl:for-each>
	    		</xsl:for-each>
	    		</tr>
                <xsl:for-each select="sample">
                    <xsl:variable name="prop_filter" select="@proportion_filter * 100"/>                
                    <tr>
                        <td nowrap="true"><a href="{/analysis/@run_name}/{./@name}.html"><xsl:value-of select="@name"/></a></td>
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
                    <xsl:apply-templates select="."/>
                </xsl:for-each>
            </table>
	    	<em>Percentages indicate the percentage of the sample containing that mutation, a value of 'none' indicates that no resistant mutations were present in that gene.</em>
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
        	<center><h1>TB Detailed ASAP Run Summary for: <xsl:value-of select="@run_name"/></h1></center>
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
	    		<th class="rotate"><div><span>Rifampin Resistance</span></div></th>
	    		<th class="rotate"><div><span>Isoniazid Resistance</span></div></th>
	    		<th class="rotate"><div><span>Quinolone Resistance</span></div></th>
	    		<th class="rotate"><div><span>Kanamycin Resistance</span></div></th>
	    		<th class="rotate"><div><span>Capreomycin Resistance</span></div></th>
	    		<th class="rotate"><div><span>Amikacin Resistance</span></div></th>
	    		</tr>
                <tr>
                    <td class="norotate"><xsl:value-of select="@name"/></td>
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
	</exsl:document>

<!-- Per Sample Detailed Results -->
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
