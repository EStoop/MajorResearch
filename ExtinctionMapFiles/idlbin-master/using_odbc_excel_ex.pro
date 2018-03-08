;
;IDL Read Excel file through ODBC
;
;Author: DongYanqing
;dongyq@esrichina-bj.cn
;http://hi.baidu.com/dyqwrp
;
PRO using_odbc_excel_ex
  filename = file_dirname(ROUTINE_FILEPATH('Using_ODBC_EXCEL'))+'\data.xlsx'
  ;Check DB
  IF DB_EXISTS() EQ 0 THEN BEGIN
    msg = DIALOG_MESSAGE('²»Ö§³ÖODBC!',/Error)
    RETURN
  ENDIF
  ;Create new DB object
  oDatabase = OBJ_NEW('IDLdbDatabase')
  ;check the DataSources
  sources = oDatabase->GETDATASOURCES()
  index = WHERE(sources.DATASOURCE EQ 'Excel Files',count)
  IF count EQ 0 THEN BEGIN
    msg = DIALOG_MESSAGE('ODBC can not read Excel Files',/Error)
    OBJ_DESTROY,oDatabase
    RETURN
  ENDIF
  ;check file exist
  IF ~FILE_TEST(filename) THEN BEGIN
    msg = DIALOG_MESSAGE('File doesnot exist!',/Error)
    OBJ_DESTROY,oDatabase
    RETURN
  ENDIF
  ;connect
  oDatabase->CONNECT,DATASOURCE='Excel Files;DBQ='+filename
  ;check the connection 
  oDatabase->GETPROPERTY,IS_CONNECTED = connectStat
  IF connectStat EQ 0 THEN BEGIN
    msg = DIALOG_MESSAGE('DB connection unSuccessful...',/Error)
    OBJ_DESTROY,oDatabase
    RETURN
  ENDIF
  
  ;Get Tables
  tables = oDatabase->GETTABLES()
  nTables = N_ELEMENTS(tables)
  FOR i=0, nTables-1 DO BEGIN
    ;Make sure that the table name should be added []
        tname = '[' + tables[i].NAME + ']' 
    PRINT, 'table name', tname
    
    oRecordset = OBJ_NEW('IDLdbRecordset',oDatabase,table=tname);, SQL=sqlstr)
    ;get the field information
    oRecordset->GETPROPERTY,field_info = fieldinfo
    NFileds = N_ELEMENTS(fieldinfo)
    ;moave and read all records
    IF oRecordset->MOVECURSOR(/first) THEN BEGIN      
      FOR j=0, NFileds-1 DO BEGIN
        Value = oRecordset->GETFIELD(j)
        PRINT, 'Talbe: ' + (fieldinfo.TABLE_NAME)[j] + ', ' + $
          'Filed Name: ' + (fieldinfo.FIELD_NAME)[j] + ', ' + $
          'Value: ', Value
      ENDFOR      
      WHILE oRecordset->MOVECURSOR(/next) DO BEGIN        
        FOR j=0, NFileds-1 DO BEGIN
          Value = oRecordset->GETFIELD(j)
          PRINT, 'Talbe: ' + (fieldinfo.TABLE_NAME)[j] + ', ' + $
            'Filed Name: ' + (fieldinfo.FIELD_NAME)[j] + ', ' + $
            'Value: ', Value
        ENDFOR
      ENDWHILE
    ENDIF ELSE BEGIN 
      msg = DIALOG_MESSAGE('No Records!',/Infor)
      OBJ_DESTROY, oRecordset      
    ENDELSE    
    ;Destroy the object
    OBJ_DESTROY, oRecordset
  ENDFOR
  OBJ_DESTROY,oDatabase
END