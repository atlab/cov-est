function obj = getSchema
persistent schemaObject

if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'sim', 'dimitri_simulation');
end

obj = schemaObject;
end