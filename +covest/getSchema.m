function schema = getSchema
persistent schemaObject

if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'covest', 'dimitri_covest');
end

schema = schemaObject;
end