Imports ESRI.ArcGIS.esriSystem
Imports System.IO

Public Class clsExtensionTT
    Inherits ESRI.ArcGIS.Desktop.AddIns.Extension

    Private _ExtensionState As esriExtensionState
    Private Shared s_extension As clsExtensionTT

    Public Sub New()
        s_extension = Me
    End Sub

    Protected Overrides Sub OnStartup()

    End Sub

    Protected Overrides Sub OnShutdown()

    End Sub

    Friend Shared Function GetExtension()
        ' Extension loads just in time, call FindExtension to load it.
        If s_extension Is Nothing Then
            Dim extID As UID = New UIDClass()
            extID.Value = My.ThisAddIn.IDs.clsExtensionTT
            My.ArcMap.Application.FindExtensionByCLSID(extID)
        End If
        Return s_extension
    End Function
End Class