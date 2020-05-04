Imports ESRI.ArcGIS.Geoprocessing
Imports ESRI.ArcGIS.GeoprocessingUI
Imports ESRI.ArcGIS.esriSystem
Imports ESRI.ArcGIS.Geodatabase

Public Class clsTravelTime
    Inherits ESRI.ArcGIS.Desktop.AddIns.Button
    Private _Settings As clsExtensionTT = Nothing

    Public Sub New()

    End Sub

    Protected Overrides Sub OnClick()
        Try
            Dim pExt As clsExtensionTT = clsExtensionTT.GetExtension()

            '***Get the settings from the extension
            If pExt Is Nothing Then
                MsgBox("The Travel Time extension was not found.", MsgBoxStyle.Exclamation, MSG_TITLE)
                Exit Sub
            End If

            'Execute Python Script
            ExecutePythonTool("traveltime")
        Catch ex As Exception
            LogError(ex, True)
        End Try
    End Sub

    Protected Overrides Sub OnUpdate()

    End Sub

    Private Sub ExecutePythonTool(ByVal strTool As String)
        Try
            'Set a reference to the IGPCommandHelper2 interface.
            Dim pToolHelper As IGPToolCommandHelper2 = New GPToolCommandHelper

            'Set the tool you want to invoke.
            Dim toolboxPath As String = vbNullString
            toolboxPath = System.IO.Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location) & "\Files\HEItraveltime.tbx"
            pToolHelper.SetToolByName(toolboxPath, strTool)

            Dim pParams As IArray
            pParams = pToolHelper.Tool.ParameterInfo

            'Invoke the tool.
            pToolHelper.Invoke(pParams)
            My.ArcMap.Application.CurrentTool = Nothing

        Catch ex As Exception
            LogError(ex, True)
        End Try
    End Sub
End Class
