function mongoTest()

%system('ssh -fN -L27017:localhost:27017 brainflight');
mongoDB = Mongo();
db = 'manuel';
collections = mongoDB.getDatabaseCollections(db);
%%
buf = BsonBuffer();
buf.append('game', 'b4b_1');
cursor = MongoCursor(buf.finish);
if mongoDB.find(collections{15}, cursor);
    idx = 1;
    while cursor.next()
        id(idx) = item.value('solution');
        idx = idx + 1;
    end
end
buf = BsonBuffer();
buf.append('game', 'b4b_1');
cursor = MongoCursor(buf.finish);

delete(cursor);

end