<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:exsl="http://exslt.org/common" xmlns:str="http://exslt.org/strings" xmlns:asap="http://pathogen.tgen.org/ASAP/functions" extension-element-prefixes="exsl str asap">
    <xsl:import href="http://exslt.org/str/functions/replace/str.replace.function.xsl"/>
    <xsl:output method="xhtml" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" omit-xml-declaration="yes" encoding="UTF-8" indent="yes"/>
    <xsl:template match="/analysis">
        <html>
        <head>
            <title>Run Summary for: <xsl:value-of select="@run_name"/></title>
            <style>
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
                }
                .second-header th{
                    top: 30px;
                }
                td {
                    border: 1px solid #ccc;
                    background-color: #fff;
                }
                .fixed-header {
                    z-index: 50;
                }
                .fixed-col-sample {
                    left: 0;
                    position: sticky;
                    width: 240px;
                    border-right: 3px solid #ccc
                }
            </style>
        </head>
        <body>
            <center><h1>ASAP Run Summary for: <xsl:value-of select="@run_name"/></h1></center>
	    <br/>
            <em>Run stats</em>
            <table>
                <tr>
                <td>Total Samples</td>
                <td><xsl:value-of select="count(//sample[not(contains(@name, 'NTC') or contains(@name, 'NEC') or contains(@name, 'PTC'))])"/></td>
                </tr>
                <tr>
                <td>Total (&gt;90%)</td>
                <td><xsl:value-of select="count(//sample[not(contains(@name, 'NTC') or contains(@name, 'NEC') or contains(@name, 'PTC'))]//breadth[text() &gt;= 90])"/></td>
                </tr>
                <tr>
                <td>Total (80-90%)</td>
                <td><xsl:value-of select="count(//sample[not(contains(@name, 'NTC') or contains(@name, 'NEC') or contains(@name, 'PTC'))]//breadth[text() &gt;= 80 and text() &lt; 90])"/></td>
                </tr>
                <tr>
                <td>Total (70-80%)</td>
                <td><xsl:value-of select="count(//sample[not(contains(@name, 'NTC') or contains(@name, 'NEC') or contains(@name, 'PTC'))]//breadth[text() &gt;= 70 and text() &lt; 80])"/></td>
                </tr>
                <tr>
                <td>Total (&lt;70%)</td>
                <td><xsl:value-of select="count(//sample[not(contains(@name, 'NTC') or contains(@name, 'NEC') or contains(@name, 'PTC'))]//breadth[text() &lt; 70])"/></td>
                </tr>
            </table>
            <br/>
	    <em>Alignment stats for each sample</em>
            <table class="freeze-table" width="100%">
                <thead>
	        <tr>
	    	<th class="fixed-col-sample fixed-header">Sample</th>
                <th>Total Reads</th>
	    	<xsl:for-each select="sample[1]">
	    	    <xsl:for-each select="assay">
	    	        <th><xsl:value-of select='@name'/> Aligned</th>
                        <th>% Aligned</th>
                        <th>Coverage Breadth</th>
                        <th>Average Depth</th>
                        <th># SNPs + INDELs (10-80%)</th>
                        <th># SNPs + INDELs (>=80%)</th>
	    	    </xsl:for-each>
	    	</xsl:for-each>
	    	</tr>
                </thead>
                <xsl:for-each select="sample">
                    <tr>
                    <td nowrap="true" class="fixed-col-sample"><a href="{/analysis/@run_name}/{./@name}.html"><xsl:value-of select="@name"/></a></td>
                    <xsl:variable name="TOTAL_READS" select="@mapped_reads + @unmapped_reads"/>
                    <td align="right"><xsl:value-of select="$TOTAL_READS"/></td>
		    <xsl:for-each select="assay">
		        <td align="right"><xsl:value-of select="amplicon/@reads"/></td>
	    		<td align="right"><xsl:value-of select='format-number(amplicon/@reads div $TOTAL_READS * 100, "##.##")'/>%</td>
	    		<td align="right"><xsl:value-of select='format-number(amplicon/breadth, "##.##")'/>%</td>
	    		<td align="right"><xsl:value-of select='format-number(amplicon/average_depth, "##.##")'/></td>
                        <td align="right"><xsl:value-of select='count(amplicon/snp/snp_call[@percent &gt;= 10 and @percent &lt; 80])'/></td>
                        <td align="right"><xsl:value-of select='count(amplicon/snp/snp_call[@percent &gt;= 80])'/></td>
		    </xsl:for-each>
                    </tr>
                    <xsl:apply-templates select="."/>
                </xsl:for-each>
            </table>
            <br/>
            <a href="{@run_name}_mutations.html">Click here for the mutation report</a><br/>
            <a href="{@run_name}_variants.html">Click here for the variant report</a><br/>
            <a href="{@run_name}_omicron.html">Click here for the omicron report</a>
        </body>
        </html>

<!-- Mutation Report -->
    <exsl:document method="html" href="{@run_name}_omicron.html">
        <html>
        <head>
            <title>Omicron Report for: <xsl:value-of select="@run_name"/></title>
            <style>
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
                }
                .second-header th{
                    top: 30px;
                }
                td {
                    border: 1px solid #ccc;
                    background-color: #fff;
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
            <center><h1>SARS-CoV-2 Omicron Mutation Report for: <xsl:value-of select="@run_name"/></h1></center>
	    <br/>
            <table class="freeze-table" width="100%">
                <thead>
	        <tr>
	    	<th rowspan="2" class="fixed-col-sample fixed-header">Sample (Total: <xsl:value-of select="count(//sample[not(contains(@name, 'NTC'))])"/>)</th>
	    	<th rowspan="2" class="fixed-col-breadth fixed-header">Coverage Breadth</th>
                <th colspan="22">All Omicron</th>
                <th colspan="11">BA.1 Only</th>
                <th colspan="6">BA.2 Only</th>
                </tr>
                <tr class="second-header">
                <th>S:G142D</th>
                <th>S:G339D</th>
                <th>S:S373P</th>
                <th>S:S375F</th>
                <th>S:K417N</th>
                <th>S:N440K</th>
                <th>S:S477N</th>
                <th>S:T478K</th>
                <th>S:E484A</th>
                <th>S:Q493R</th>
                <th>S:Q498R</th>
                <th>S:N501Y</th>
                <th>S:Y505H</th>
                <th>S:D614G</th>
                <th>S:H655Y</th>
                <th>S:N679K</th>
                <th>S:P681H</th>
                <th>S:N764K</th>
                <th>S:D796Y</th>
                <th>S:Q954H</th>
                <th>S:N969K</th>
                <th>M:Q19E</th>
                <th>S:A67V</th>
                <th>S:T95I</th>
                <th>S:N211I</th>
                <th>S:215EPEins</th>
                <th>S:S371L</th>
                <th>S:G446S</th>
                <th>S:G496S</th>
                <th>S:T547K</th>
                <th>S:N856K</th>
                <th>S:L981F</th>
                <th>M:D3G</th>
                <th>S:T19I</th>
                <th>S:V213G</th>
                <th>S:S371F</th>
                <th>S:T376A</th>
                <th>S:D405N</th>
                <th>S:R408S</th>
                </tr>
                </thead>
                <xsl:for-each select="sample">
                    <xsl:variable name="TOTAL_READS" select="@mapped_reads + @unmapped_reads"/>
                    <xsl:variable name="BG_A67V"><xsl:choose><xsl:when test="descendant::snp[@name='S:A67V']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:A67V']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:A67V']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:A67V']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_T95I"><xsl:choose><xsl:when test="descendant::snp[@name='S:T95I']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T95I']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:T95I']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T95I']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_G142D"><xsl:choose><xsl:when test="descendant::snp[@name='S:G142D']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:G142D']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:G142D']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:G142D']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_N211I"><xsl:choose><xsl:when test="descendant::snp[@name='S:N211I']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N211I']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:N211I']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N211I']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_G339D"><xsl:choose><xsl:when test="descendant::snp[@name='S:G339D']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:G339D']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:G339D']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:G339D']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_S371L"><xsl:choose><xsl:when test="descendant::snp[@name='S:S371L']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:S371L']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:S371L']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:S371L']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_S373P"><xsl:choose><xsl:when test="descendant::snp[@name='S:S373P']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:S373P']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:S373P']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:S373P']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_S375F"><xsl:choose><xsl:when test="descendant::snp[@name='S:S375F']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:S375F']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:S375F']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:S375F']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_K417N"><xsl:choose><xsl:when test="descendant::snp[@name='S:K417N']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:K417N']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:K417N']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:K417N']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_N440K"><xsl:choose><xsl:when test="descendant::snp[@name='S:N440K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N440K']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:N440K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N440K']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_G446S"><xsl:choose><xsl:when test="descendant::snp[@name='S:G446S']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:G446S']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:G446S']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:G446S']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_S477N"><xsl:choose><xsl:when test="descendant::snp[@name='S:S477N']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:S477N']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:S477N']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:S477N']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_T478K"><xsl:choose><xsl:when test="descendant::snp[@name='S:T478K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T478K']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:T478K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T478K']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_E484A"><xsl:choose><xsl:when test="descendant::snp[@name='S:E484A']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:E484A']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:E484A']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:E484A']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_Q493R"><xsl:choose><xsl:when test="descendant::snp[@name='S:Q493R']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:Q493R']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:Q493R']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:Q493R']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_G496S"><xsl:choose><xsl:when test="descendant::snp[@name='S:G496S']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:G496S']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:G496S']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:G496S']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_Q498R"><xsl:choose><xsl:when test="descendant::snp[@name='S:Q498R']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:Q498R']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:Q498R']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:Q498R']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_N501Y"><xsl:choose><xsl:when test="descendant::snp[@name='S:N501Y']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N501Y']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:N501Y']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N501Y']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_Y505H"><xsl:choose><xsl:when test="descendant::snp[@name='S:Y505H']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:Y505H']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:Y505H']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:Y505H']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_T547K"><xsl:choose><xsl:when test="descendant::snp[@name='S:T547K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T547K']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:T547K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T547K']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_D614G"><xsl:choose><xsl:when test="descendant::snp[@name='S:D614G']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:D614G']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:D614G']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:D614G']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_H655Y"><xsl:choose><xsl:when test="descendant::snp[@name='S:H655Y']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:H655Y']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:H655Y']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:H655Y']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_N679K"><xsl:choose><xsl:when test="descendant::snp[@name='S:N679K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N679K']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:N679K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N679K']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_P681H"><xsl:choose><xsl:when test="descendant::snp[@name='S:P681H']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:P681H']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:P681H']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:P681H']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_N764K"><xsl:choose><xsl:when test="descendant::snp[@name='S:N764K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N764K']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:N764K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N764K']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_D796Y"><xsl:choose><xsl:when test="descendant::snp[@name='S:D796Y']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:D796Y']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:D796Y']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:D796Y']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_N856K"><xsl:choose><xsl:when test="descendant::snp[@name='S:N856K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N856K']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:N856K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N856K']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_Q954H"><xsl:choose><xsl:when test="descendant::snp[@name='S:Q954H']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:Q954H']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:Q954H']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:Q954H']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_N969K"><xsl:choose><xsl:when test="descendant::snp[@name='S:N969K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N969K']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:N969K']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:N969K']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_L981F"><xsl:choose><xsl:when test="descendant::snp[@name='S:L981F']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:L981F']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:L981F']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:L981F']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_D3G"><xsl:choose><xsl:when test="descendant::snp[@name='M:D3G']/snp_call/@percent &gt; 80 and descendant::snp[@name='M:D3G']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='M:D3G']/snp_call/@percent &gt; 80 and descendant::snp[@name='M:D3G']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_Q19E"><xsl:choose><xsl:when test="descendant::snp[@name='M:Q19E']/snp_call/@percent &gt; 80 and descendant::snp[@name='M:Q19E']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='M:Q19E']/snp_call/@percent &gt; 80 and descendant::snp[@name='M:Q19E']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_215EPEins"><xsl:choose><xsl:when test="descendant::snp[@name='S:215EPEins']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:215EPEins']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:215EPEins']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:215EPEins']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_T19I"><xsl:choose><xsl:when test="descendant::snp[@name='S:T19I']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T19I']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:T19I']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T19I']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_V213G"><xsl:choose><xsl:when test="descendant::snp[@name='S:V213G']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:V213G']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:V213G']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:V213G']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_S371F"><xsl:choose><xsl:when test="descendant::snp[@name='S:S371F']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:S371F']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:S371F']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:S371F']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_T376A"><xsl:choose><xsl:when test="descendant::snp[@name='S:T376A']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T376A']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:T376A']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:T376A']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_D405N"><xsl:choose><xsl:when test="descendant::snp[@name='S:D405N']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:D405N']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:D405N']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:D405N']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="BG_R408S"><xsl:choose><xsl:when test="descendant::snp[@name='S:R408S']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:R408S']/@depth &gt;= 50">#4EB5C8</xsl:when><xsl:when test="descendant::snp[@name='S:R408S']/snp_call/@percent &gt; 80 and descendant::snp[@name='S:R408S']/@depth &gt;= 10">#B1DBE3</xsl:when><xsl:otherwise>#FFFFFF</xsl:otherwise></xsl:choose></xsl:variable>
                    <xsl:variable name="COLOR"><xsl:choose><xsl:when test="assay/amplicon/breadth &lt; 90">#888888</xsl:when><xsl:otherwise>#000000</xsl:otherwise></xsl:choose></xsl:variable>
                    <tr>
                    <td nowrap="true" class="fixed-col-sample"><a href="{/analysis/@run_name}/{./@name}.html"><xsl:value-of select="@name"/></a></td>
	    	    <td align="center" style="color:{$COLOR}" class="fixed-col-breadth"><xsl:value-of select='format-number(assay/amplicon/breadth, "##.##")'/>%</td>
                    <td style="background:{$BG_G142D}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:G142D']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:G142D']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:G142D"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_G339D}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:G339D']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:G339D']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:G339D"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_S373P}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:S373P']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:S373P']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:S373P"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_S375F}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:S375F']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:S375F']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:S375F"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_K417N}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:K417N']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:K417N']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:K417N"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_N440K}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:N440K']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:N440K']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:N440K"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_S477N}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:S477N']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:S477N']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:S477N"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_T478K}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:T478K']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:T478K']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:T478K"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_E484A}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:E484A']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:E484A']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:E484A"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_Q493R}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:Q493R']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:Q493R']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:Q493R"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_Q498R}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:Q498R']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:Q498R']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:Q498R"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_N501Y}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:N501Y']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:N501Y']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:N501Y"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_Y505H}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:Y505H']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:Y505H']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:Y505H"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_D614G}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:D614G']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:D614G']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:D614G"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_H655Y}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:H655Y']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:H655Y']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:H655Y"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_N679K}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:N679K']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:N679K']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:N679K"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_P681H}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:P681H']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:P681H']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:P681H"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_N764K}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:N764K']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:N764K']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:N764K"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_D796Y}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:D796Y']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:D796Y']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:D796Y"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_Q954H}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:Q954H']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:Q954H']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:Q954H"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_N969K}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:N969K']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:N969K']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:N969K"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_Q19E}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='M:Q19E']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='M:Q19E']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="M:Q19E"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_A67V}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:A67V']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:A67V']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:A67V"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_T95I}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:T95I']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:T95I']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:T95I"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_N211I}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:N211I']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:N211I']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:N211I"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_215EPEins}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:215EPEins']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:215EPEins']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:215EPEins"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_S371L}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:S371L']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:S371L']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:S371L"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_G446S}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:G446S']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:G446S']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:G446S"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_G496S}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:G496S']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:G496S']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:G496S"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_T547K}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:T547K']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:T547K']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:T547K"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_N856K}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:N856K']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:N856K']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:N856K"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_L981F}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:L981F']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:L981F']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:L981F"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_D3G}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='M:D3G']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='M:D3G']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="M:D3G"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_T19I}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:T19I']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:T19I']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:T19I"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_V213G}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:V213G']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:V213G']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:V213G"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_S371F}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:S371F']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:S371F']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:S371F"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_T376A}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:T376A']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:T376A']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:T376A"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_D405N}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:D405N']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:D405N']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:D405N"]/snp_call/@percent, "##.##")'/>%)</td>
                    <td style="background:{$BG_R408S}; color:{$COLOR}"><xsl:value-of select="descendant::snp[@name='S:R408S']/snp_call/@count"/>/<xsl:value-of select="descendant::snp[@name='S:R408S']/@depth"/>(<xsl:value-of select='format-number(descendant::snp[@name="S:R408S"]/snp_call/@percent, "##.##")'/>%)</td>
                    </tr>
                </xsl:for-each>
            </table>
            <br/>
            <a href="{@run_name}.html">Click here to return to the run summary</a>
        </body>
        </html>

    </exsl:document>
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
	        Total reads: <xsl:value-of select="@mapped_reads + @unmapped_reads"/><br/>
	        Mapped reads: <xsl:value-of select="@mapped_reads"/><br/>
	        Unmapped reads: <xsl:value-of select="@unmapped_reads"/><br/>
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
	    	<h3>Mutations present in sample: <xsl:value-of select="@name"/></h3>
	    	<table border="2" width="100%">
	    		<tr>
	    		<th>Assay Name</th>
	    		<th># of Reads</th>
	    		<th>Coverage Breadth</th>
	    		<th>Average Depth</th>
	    		<th>Tracked Mutations</th>
	    		<th>Other SNPs found(% reads containing SNP)</th>
	    		</tr>
	    		<xsl:for-each select="assay">
	    		    <xsl:if test="@type = 'SNP' or @type = 'mixed' and amplicon/@reads &gt; 0">
	    		    <tr>
	    		        <td><xsl:value-of select="@name"/></td>
	    		        <td><xsl:value-of select="amplicon/@reads"/></td>
	    		        <td><xsl:value-of select='format-number(amplicon/breadth, "##.##")'/>%</td>
	    		        <td><xsl:value-of select='format-number(amplicon/average_depth, "##.###")'/></td>
	    		        <td><xsl:for-each select="amplicon/snp">
                                    <xsl:if test="@name!='unknown'">
	    		            <xsl:value-of select="@name"/> - <xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%)<br/>
                                    </xsl:if>
	    		        </xsl:for-each></td>
	    		        <!--<td><xsl:value-of select="amplicon/significance"/><xsl:if test="amplicon/significance/@flag"> (<xsl:value-of select="amplicon/significance/@flag"/>)</xsl:if></td>-->
	    		        <td><xsl:for-each select="amplicon/snp">
                                    <xsl:if test="@name='unknown'">
	    		            <xsl:value-of select="@position"/><xsl:value-of select="@reference"/>-><xsl:value-of select="snp_call"/>(<xsl:value-of select='format-number(snp_call/@percent, "##.##")'/>%) <xsl:value-of select="snp_call/@count"/>/<xsl:value-of select="@depth"/><br/>
                                    </xsl:if>
	    		        </xsl:for-each></td>
	    		    </tr>
	    		    </xsl:if>
	    		</xsl:for-each>
	    	</table>
	    	<br />
	    </body>
	    </html>
	</exsl:document>
<xsl:variable name="SAMPLE">
<xsl:choose><xsl:when test="contains(@name, 'WMTS')">
<xsl:value-of select="substring(@name,1,25)"/>
</xsl:when><xsl:when test="contains(@name, 'ARTIC')">
<xsl:value-of select="substring(@name,1,26)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled-NEC') or contains(@name, 'Tiled_NEC')">
<xsl:value-of select="substring(@name,1,33)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled-TG1') or contains(@name, 'Tiled_TG1')">
<xsl:value-of select="substring(@name,1,27)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled-') or contains(@name, 'Tiled_')">
<xsl:value-of select="substring(@name,1,26)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled400')">
<xsl:value-of select="substring(@name,1,29)"/>
</xsl:when><xsl:when test="contains(@name, 'Tiled1000')">
<xsl:value-of select="substring(@name,1,30)"/>
</xsl:when><xsl:otherwise>
<xsl:value-of select="@name"/>
</xsl:otherwise></xsl:choose>
</xsl:variable>
  <xsl:if test="assay[1]/amplicon[1]/breadth &gt;= 90 and assay[1]/amplicon[1]/average_depth &gt;= 30">
<exsl:document method="text" href="{/analysis/@run_name}/{$SAMPLE}_gapfilled.fasta">
<xsl:for-each select="assay">
<xsl:text/>><xsl:value-of select="$SAMPLE"/>
<xsl:text>&#xa;</xsl:text>
<xsl:value-of select="./amplicon[1]/gapfilled_consensus_sequence"/>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
<xsl:text/>
</exsl:document>
  </xsl:if>
  <xsl:if test="assay[1]/amplicon[1]/breadth &gt;= 80 and assay[1]/amplicon[1]/breadth &lt; 90 and assay[1]/amplicon[1]/average_depth &gt;= 30">
<exsl:document method="text" href="{/analysis/@run_name}/{$SAMPLE}_incomplete.fasta">
<xsl:for-each select="assay">
<xsl:text/>><xsl:value-of select="$SAMPLE"/>
<xsl:text>&#xa;</xsl:text>
<xsl:value-of select="./amplicon[1]/gapfilled_consensus_sequence"/>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
<xsl:text/>
</exsl:document>
  </xsl:if>
  <xsl:if test="assay[1]/amplicon[1]/breadth &lt; 80 or assay[1]/amplicon[1]/average_depth &lt; 30">
<xsl:variable name="BREADTH" select='format-number(assay/amplicon/breadth, "0")'/>
<exsl:document method="text" href="{/analysis/@run_name}/{$SAMPLE}_{$BREADTH}percent.fasta">
<xsl:for-each select="assay">
<xsl:text/>><xsl:value-of select="$SAMPLE"/>
<xsl:text>&#xa;</xsl:text>
<xsl:value-of select="./amplicon[1]/gapfilled_consensus_sequence"/>
<xsl:text>&#xa;</xsl:text>
</xsl:for-each>
<xsl:text/>
</exsl:document>
  </xsl:if>
    </xsl:template>
</xsl:stylesheet>
