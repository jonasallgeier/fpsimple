    (function() {
          var fn = function() {
            Bokeh.safely(function() {
              (function(root) {
                function embed_document(root) {
                  
                var docs_json = '{"aad0f763-9f6d-4b1d-a19f-82e554a0c78f":{"defs":[],"roots":{"references":[{"attributes":{"bar_color":"rgba(60, 60, 60, 0.6)","end":5.0,"format":{"id":"1143"},"js_property_callbacks":{"change:value":[{"id":"1160"}]},"sizing_mode":"stretch_width","start":0.1,"step":0.1,"value":1.0},"id":"1144","type":"Slider"},{"attributes":{"data":{"x":{"__ndarray__":"AAAAAAAAAAAa60NjfWg0QBrrQ2N9aERAp+DlFLycTkAa60NjfWhUQODlFLycgllAp+DlFLycXkC3bdu2bdthQBrrQ2N9aGRAfWisD431ZkDg5RS8nIJpQERjfWisD2xAp+DlFLycbkAFL6fg5ZRwQLdt27Zt23FAaKwPjfUhc0Aa60NjfWh0QMwpeDkFr3VAfWisD431dkAvp+DlFDx4QODlFLycgnlAkiRJkiTJekBEY31orA98QPWhsT40Vn1Ap+DlFLycfkBZHxrrQ+N/QAUvp+DllIBAXk7Byyk4gUC3bdu2bduBQBCN9aGxfoJAaKwPjfUhg0DByyl4OcWDQBrrQ2N9aIRAcwpeTsELhUDMKXg5Ba+FQCRJkiRJUoZAfWisD431hkDWh8b60JiHQC+n4OUUPIhAiMb60FjfiEDg5RS8nIKJQDkFL6fgJYpAkiRJkiTJikDrQ2N9aGyLQERjfWisD4xAnYKXU/CyjED1obE+NFaNQE7Byyl4+Y1Ap+DlFLycjkAAAAAAAECPQA==","dtype":"float64","order":"little","shape":[50]},"yBot":{"__ndarray__":"AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fw==","dtype":"float64","order":"little","shape":[50]},"yTop":{"__ndarray__":"AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fw==","dtype":"float64","order":"little","shape":[50]}},"selected":{"id":"1179"},"selection_policy":{"id":"1180"}},"id":"1078","type":"ColumnDataSource"},{"attributes":{"data_source":{"id":"1078"},"glyph":{"id":"1095"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"1096"},"view":{"id":"1098"}},"id":"1097","type":"GlyphRenderer"},{"attributes":{"fill_alpha":0.4,"fill_color":"rgb(105, 186, 201)","x":{"field":"x"},"y1":{"field":"yTop"},"y2":{"field":"yBot"}},"id":"1095","type":"VArea"},{"attributes":{"sizing_mode":"stretch_width","style":{"font-size":"115%"},"text":"&lt;i&gt;L&lt;/i&gt;"},"id":"1145","type":"Div"},{"attributes":{"data":{"x":{"__ndarray__":"AAAAAAAAAAAa60NjfWg0QBrrQ2N9aERAp+DlFLycTkAa60NjfWhUQODlFLycgllAp+DlFLycXkC3bdu2bdthQBrrQ2N9aGRAfWisD431ZkDg5RS8nIJpQERjfWisD2xAp+DlFLycbkAFL6fg5ZRwQLdt27Zt23FAaKwPjfUhc0Aa60NjfWh0QMwpeDkFr3VAfWisD431dkAvp+DlFDx4QODlFLycgnlAkiRJkiTJekBEY31orA98QPWhsT40Vn1Ap+DlFLycfkBZHxrrQ+N/QAUvp+DllIBAXk7Byyk4gUC3bdu2bduBQBCN9aGxfoJAaKwPjfUhg0DByyl4OcWDQBrrQ2N9aIRAcwpeTsELhUDMKXg5Ba+FQCRJkiRJUoZAfWisD431hkDWh8b60JiHQC+n4OUUPIhAiMb60FjfiEDg5RS8nIKJQDkFL6fgJYpAkiRJkiTJikDrQ2N9aGyLQERjfWisD4xAnYKXU/CyjED1obE+NFaNQE7Byyl4+Y1Ap+DlFLycjkAAAAAAAECPQA==","dtype":"float64","order":"little","shape":[50]},"yBot":{"__ndarray__":"AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fw==","dtype":"float64","order":"little","shape":[50]},"yTop":{"__ndarray__":"AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fwAAAAAAAPh/AAAAAAAA+H8AAAAAAAD4fw==","dtype":"float64","order":"little","shape":[50]}},"selected":{"id":"1181"},"selection_policy":{"id":"1182"}},"id":"1079","type":"ColumnDataSource"},{"attributes":{"data":{"x":{"__ndarray__":"AAAAAAAAAAAa60NjfWg0QBrrQ2N9aERAp+DlFLycTkAa60NjfWhUQODlFLycgllAp+DlFLycXkC3bdu2bdthQBrrQ2N9aGRAfWisD431ZkDg5RS8nIJpQERjfWisD2xAp+DlFLycbkAFL6fg5ZRwQLdt27Zt23FAaKwPjfUhc0Aa60NjfWh0QMwpeDkFr3VAfWisD431dkAvp+DlFDx4QODlFLycgnlAkiRJkiTJekBEY31orA98QPWhsT40Vn1Ap+DlFLycfkBZHxrrQ+N/QAUvp+DllIBAXk7Byyk4gUC3bdu2bduBQBCN9aGxfoJAaKwPjfUhg0DByyl4OcWDQBrrQ2N9aIRAcwpeTsELhUDMKXg5Ba+FQCRJkiRJUoZAfWisD431hkDWh8b60JiHQC+n4OUUPIhAiMb60FjfiEDg5RS8nIKJQDkFL6fgJYpAkiRJkiTJikDrQ2N9aGyLQERjfWisD4xAnYKXU/CyjED1obE+NFaNQE7Byyl4+Y1Ap+DlFLycjkAAAAAAAECPQA==","dtype":"float64","order":"little","shape":[50]},"y":{"__ndarray__":"AAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQAAAAAAAwHJAAAAAAADAckAAAAAAAMByQA==","dtype":"float64","order":"little","shape":[50]}},"selected":{"id":"1173"},"selection_policy":{"id":"1174"}},"id":"1077","type":"ColumnDataSource"},{"attributes":{"below":[{"id":"1011"}],"center":[{"id":"1014"},{"id":"1018"}],"height_policy":"fit","left":[{"id":"1015"}],"margin":[5,0,0,0],"match_aspect":true,"min_width":250,"renderers":[{"id":"1087"},{"id":"1092"},{"id":"1097"},{"id":"1102"},{"id":"1107"},{"id":"1112"},{"id":"1117"},{"id":"1122"}],"title":{"id":"1055"},"toolbar":{"id":"1023"},"width_policy":"max","x_range":{"id":"1003"},"x_scale":{"id":"1007"},"y_range":{"id":"1005"},"y_scale":{"id":"1009"}},"id":"1002","subtype":"Figure","type":"Plot"},{"attributes":{"data":{"f":[0,"NaN",100],"t":[0,50,100]},"selected":{"id":"1185"},"selection_policy":{"id":"1186"}},"id":"1080","type":"ColumnDataSource"},{"attributes":{},"id":"1009","type":"LinearScale"},{"attributes":{"data":{"clr":[[0.9],[0.8],[0.7],[0.6],[0.5],[0.4],[0.30000000000000004],[0.19999999999999996],[0.09999999999999998]],"xs":[[100.0,100.0],[200.0,200.0],[300.0,300.0],[400.0,400.0],[500.0,500.0],[600.0,600.0],[700.0,700.0],[800.0,800.0],[900.0,900.0]],"ys":[[0,300.0],[0,300.0],[0,300.0],[0,300.0],[0,300.0],[0,300.0],[0,300.0],[0,300.0],[0,300.0]]},"selected":{"id":"1177"},"selection_policy":{"id":"1178"}},"id":"1081","type":"ColumnDataSource"},{"attributes":{},"id":"1003","type":"DataRange1d"},{"attributes":{"axis_label":"x","axis_label_text_font_size":"25px","formatter":{"id":"1062"},"major_label_policy":{"id":"1061"},"major_label_text_font_size":"25px","ticker":{"id":"1012"}},"id":"1011","type":"LinearAxis"},{"attributes":{"data":{"xs":[[0,1000],[0,1000],[0,1000],[0,1000],[0,1000],[0,1000],[0,1000],[0,1000],[0,1000]],"ys":[[30.0,30.0],[60.0,60.0],[90.0,90.0],[120.0,120.0],[150.0,150.0],[180.0,180.0],[210.0,210.0],[240.0,240.0],[270.0,270.0]]},"selected":{"id":"1183"},"selection_policy":{"id":"1184"}},"id":"1082","type":"ColumnDataSource"},{"attributes":{"line_alpha":{"value":0.1},"line_color":{"field":"clr","transform":{"id":"1083"}},"line_width":{"value":5},"xs":{"field":"xs"},"ys":{"field":"ys"}},"id":"1101","type":"MultiLine"},{"attributes":{"high":1,"low":0,"palette":["#2050da","#2551db","#2952db","#2d53db","#3155db","#3456db","#3857dc","#3b58dc","#3d59dc","#405adc","#435bdd","#455cdd","#485ddd","#4a5fdd","#4d60dd","#4f61de","#5162de","#5363de","#5664de","#5865de","#5a67df","#5c68df","#5e69df","#606adf","#616be0","#636ce0","#656ee0","#676fe0","#6970e0","#6b71e1","#6c72e1","#6e74e1","#7075e1","#7176e1","#7377e2","#7578e2","#7679e2","#787be2","#7a7ce2","#7b7de2","#7d7ee3","#7e7fe3","#8081e3","#8282e3","#8383e3","#8584e4","#8686e4","#8887e4","#8988e4","#8b89e4","#8c8ae5","#8e8ce5","#8f8de5","#908ee5","#928fe5","#9391e5","#9592e6","#9693e6","#9894e6","#9996e6","#9a97e6","#9c98e7","#9d99e7","#9e9be7","#a09ce7","#a19de7","#a39ee7","#a4a0e8","#a5a1e8","#a7a2e8","#a8a3e8","#a9a5e8","#aba6e8","#aca7e9","#ada8e9","#afaae9","#b0abe9","#b1ace9","#b2aee9","#b4afea","#b5b0ea","#b6b1ea","#b8b3ea","#b9b4ea","#bab5ea","#bbb7ea","#bdb8eb","#beb9eb","#bfbaeb","#c0bceb","#c2bdeb","#c3beeb","#c4c0ec","#c5c1ec","#c7c2ec","#c8c4ec","#c9c5ec","#cac6ec","#ccc8ec","#cdc9ec","#cecaed","#cfcbed","#d1cded","#d2ceed","#d3cfed","#d4d1ed","#d5d2ed","#d7d3ee","#d8d5ee","#d9d6ee","#dad7ee","#dbd9ee","#dddaee","#dedbee","#dfddee","#e0deee","#e1dfee","#e3e0ee","#e4e1ee","#e5e2ee","#e6e3ee","#e7e4ee","#e8e5ed","#e9e6ed","#eae6ec","#ebe7ec","#ece7eb","#ece7ea","#ede7e9","#eee6e8","#eee6e6","#efe5e5","#efe5e4","#f0e4e2","#f0e3e0","#f0e2df","#f1e0dd","#f1dfdb","#f1ded9","#f1dcd7","#f1dbd6","#f1d9d4","#f1d8d2","#f1d6d0","#f1d5ce","#f1d3cc","#f1d2ca","#f1d0c8","#f1cfc6","#f1cdc4","#f1ccc2","#f1cac0","#f1c9bf","#f1c7bd","#f1c5bb","#f1c4b9","#f0c2b7","#f0c1b5","#f0bfb3","#f0beb1","#f0bcaf","#f0bbad","#efb9ac","#efb7aa","#efb6a8","#efb4a6","#efb3a4","#eeb1a2","#eeb0a0","#eeae9e","#eead9d","#edab9b","#eda999","#eda897","#eda695","#eca593","#eca391","#eca290","#eba08e","#eb9f8c","#eb9d8a","#ea9b88","#ea9a87","#ea9885","#e99783","#e99581","#e9947f","#e8927d","#e8907c","#e88f7a","#e78d78","#e78c76","#e68a75","#e68973","#e58771","#e5856f","#e5846d","#e4826c","#e4816a","#e37f68","#e37d66","#e27c65","#e27a63","#e17961","#e1775f","#e0755e","#e0745c","#df725a","#df7158","#de6f57","#de6d55","#dd6c53","#dd6a52","#dc6850","#db674e","#db654c","#da634b","#da6249","#d96047","#d95e46","#d85d44","#d75b42","#d75940","#d6573f","#d5563d","#d5543b","#d4523a","#d45038","#d34e36","#d24d35","#d24b33","#d14931","#d04730","#d0452e","#cf432c","#ce412a","#ce3f29","#cd3d27","#cc3b25","#cb3924","#cb3722","#ca3420","#c9321e","#c9301d","#c82d1b","#c72b19","#c62817","#c62515","#c52214","#c41f12","#c31c10","#c3180e","#c2140c","#c10f09","#c00907","#bf0205"]},"id":"1083","type":"LinearColorMapper"},{"attributes":{"line_color":{"field":"clr","transform":{"id":"1083"}},"line_width":{"value":5},"xs":{"field":"xs"},"ys":{"field":"ys"}},"id":"1100","type":"MultiLine"},{"attributes":{"below":[{"id":"1037"}],"center":[{"id":"1040"},{"id":"1044"}],"height_policy":"fit","left":[{"id":"1041"}],"margin":[5,0,0,0],"min_width":250,"renderers":[{"id":"1127"}],"title":{"id":"1203"},"toolbar":{"id":"1049"},"width_policy":"max","x_range":{"id":"1029"},"x_scale":{"id":"1033"},"y_range":{"id":"1031"},"y_scale":{"id":"1035"}},"id":"1028","subtype":"Figure","type":"Plot"},{"attributes":{"source":{"id":"1078"}},"id":"1098","type":"CDSView"},{"attributes":{"fill_alpha":0.1,"fill_color":"rgb(105, 186, 201)","x":{"field":"x"},"y1":{"field":"yTop"},"y2":{"field":"yBot"}},"id":"1096","type":"VArea"},{"attributes":{},"id":"1055","type":"Title"},{"attributes":{"data_source":{"id":"1076"},"glyph":{"id":"1110"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"1111"},"view":{"id":"1113"}},"id":"1112","type":"GlyphRenderer"},{"attributes":{"line_alpha":0.1,"line_cap":"round","line_color":"rgba(192,2,6,1.0)","line_width":5,"x":{"field":"xW"},"y":{"field":"yW"}},"id":"1106","type":"Line"},{"attributes":{"data_source":{"id":"1081"},"glyph":{"id":"1100"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"1101"},"view":{"id":"1103"}},"id":"1102","type":"GlyphRenderer"},{"attributes":{"source":{"id":"1081"}},"id":"1103","type":"CDSView"},{"attributes":{"args":{"hE":{"id":"1112"},"hIsoH":{"id":"1102"},"hIsoP":{"id":"1087"},"hN":{"id":"1122"},"hS":{"id":"1117"},"hW":{"id":"1107"}},"code":"\\n    var act         = this.active\\n    hW.visible      = act.includes(0)\\n    hS.visible      = act.includes(1)\\n    hE.visible      = act.includes(2)\\n    hN.visible      = act.includes(3)\\n"},"id":"1155","type":"CustomJS"},{"attributes":{"line_cap":"round","line_color":"rgba(192,2,6,1.0)","line_width":5,"x":{"field":"xW"},"y":{"field":"yW"}},"id":"1105","type":"Line"},{"attributes":{"data":{"xE":[1000,1000],"xS":[0,1000],"xW":[0,0],"yE":{"__ndarray__":"AAAAAAAAAAAAAAAAAMByQA==","dtype":"float64","order":"little","shape":[2]},"yS":[0,0],"yW":{"__ndarray__":"AAAAAAAAAAAAAAAAAMByQA==","dtype":"float64","order":"little","shape":[2]}},"selected":{"id":"1175"},"selection_policy":{"id":"1176"}},"id":"1076","type":"ColumnDataSource"},{"attributes":{"data_source":{"id":"1076"},"glyph":{"id":"1105"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"1106"},"view":{"id":"1108"}},"id":"1107","type":"GlyphRenderer"},{"attributes":{"line_cap":"round","line_color":"rgba(32, 81, 219, 1.0)","line_width":5,"x":{"field":"xE"},"y":{"field":"yE"}},"id":"1110","type":"Line"},{"attributes":{"source":{"id":"1076"}},"id":"1108","type":"CDSView"},{"attributes":{"line_cap":"round","line_color":"rgba(105,186,201,1.0)","line_width":5,"x":{"field":"xS"},"y":{"field":"yS"}},"id":"1115","type":"Line"},{"attributes":{"line_alpha":0.1,"line_cap":"round","line_color":"rgba(32, 81, 219, 1.0)","line_width":5,"x":{"field":"xE"},"y":{"field":"yE"}},"id":"1111","type":"Line"},{"attributes":{"data_source":{"id":"1076"},"glyph":{"id":"1115"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"1116"},"view":{"id":"1118"}},"id":"1117","type":"GlyphRenderer"},{"attributes":{"source":{"id":"1076"}},"id":"1113","type":"CDSView"},{"attributes":{"line_cap":"round","line_color":"rgba(201,174,105,1.0)","line_width":5,"x":{"field":"x"},"y":{"field":"y"}},"id":"1120","type":"Line"},{"attributes":{"active":[0,1,2,3],"js_property_callbacks":{"change:active":[{"id":"1155"}]},"labels":["west","south","east","north"],"margin":[10,3,10,3],"sizing_mode":"stretch_width"},"id":"1154","type":"CheckboxButtonGroup"},{"attributes":{"line_alpha":0.1,"line_cap":"round","line_color":"rgba(105,186,201,1.0)","line_width":5,"x":{"field":"xS"},"y":{"field":"yS"}},"id":"1116","type":"Line"},{"attributes":{"active":0,"js_property_callbacks":{"change:active":[{"id":"1159"},{"id":"1160"}]},"labels":["cosinusoidal","bump","composite"],"margin":[10,3,10,3],"sizing_mode":"stretch_width"},"id":"1153","type":"RadioButtonGroup"},{"attributes":{"data_source":{"id":"1077"},"glyph":{"id":"1120"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"1121"},"view":{"id":"1123"}},"id":"1122","type":"GlyphRenderer"},{"attributes":{"source":{"id":"1076"}},"id":"1118","type":"CDSView"},{"attributes":{"sizing_mode":"stretch_width","style":{"font-size":"115%"},"text":"&lt;i&gt;w&lt;sub&gt;max&lt;/sub&gt;&lt;/i&gt;/&lt;i&gt;L&lt;/i&gt;"},"id":"1151","type":"Div"},{"attributes":{"sizing_mode":"stretch_width","style":{"font-size":"115%"},"text":"&lt;i&gt;\\u03a6&lt;/i&gt;"},"id":"1152","type":"Div"},{"attributes":{"line_alpha":0.1,"line_cap":"round","line_color":"rgba(201,174,105,1.0)","line_width":5,"x":{"field":"x"},"y":{"field":"y"}},"id":"1121","type":"Line"},{"attributes":{"sizing_mode":"stretch_width","style":{"font-size":"115%"},"text":"&lt;i&gt;Q&lt;/i&gt;&lt;sup&gt;~&lt;/sup&gt;&lt;sub&gt;north&lt;/sub&gt;"},"id":"1150","type":"Div"},{"attributes":{},"id":"1058","type":"AllLabels"},{"attributes":{"data_source":{"id":"1080"},"glyph":{"id":"1125"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"1126"},"view":{"id":"1128"}},"id":"1127","type":"GlyphRenderer"},{"attributes":{"source":{"id":"1077"}},"id":"1123","type":"CDSView"},{"attributes":{},"id":"1059","type":"BasicTickFormatter"},{"attributes":{"line_cap":"round","line_color":"rgba(219,33,82,1.0)","line_width":5,"x":{"field":"t"},"y":{"field":"f"}},"id":"1125","type":"Line"},{"attributes":{"line_alpha":0.1,"line_cap":"round","line_color":"rgba(219,33,82,1.0)","line_width":5,"x":{"field":"t"},"y":{"field":"f"}},"id":"1126","type":"Line"},{"attributes":{"bar_color":"rgba(60, 60, 60, 0.6)","end":-3.506557897319982,"format":{"id":"1131"},"js_property_callbacks":{"change:value":[{"id":"1160"}]},"sizing_mode":"stretch_width","start":-8.111728083308073,"step":0.04605170185988091,"value":-5.809142990314028},"id":"1132","type":"Slider"},{"attributes":{"code":"return tick.toFixed(0)+&#x27; [L]&#x27;"},"id":"1129","type":"FuncTickFormatter"},{"attributes":{"source":{"id":"1080"}},"id":"1128","type":"CDSView"},{"attributes":{"code":"return Math.exp(tick).toExponential(1).toString()+&#x27; [L/L]&#x27;"},"id":"1131","type":"FuncTickFormatter"},{"attributes":{"bar_color":"rgba(60, 60, 60, 0.6)","end":0.5,"format":{"id":"1133"},"js_property_callbacks":{"change:value":[{"id":"1159"},{"id":"1160"}]},"sizing_mode":"stretch_width","start":0.1,"step":0.002,"value":0.3},"id":"1134","type":"Slider"},{"attributes":{"bar_color":"rgba(60, 60, 60, 0.6)","end":3000,"format":{"id":"1129"},"js_property_callbacks":{"change:value":[{"id":"1159"},{"id":"1160"}]},"sizing_mode":"stretch_width","start":100,"step":10,"value":1000},"id":"1130","type":"Slider"},{"attributes":{"bar_color":"rgba(60, 60, 60, 0.6)","end":1.0,"format":{"id":"1135"},"js_property_callbacks":{"change:value":[{"id":"1159"},{"id":"1160"}]},"sizing_mode":"stretch_width","start":0.4,"step":0.006,"value":1.0},"id":"1136","type":"Slider"},{"attributes":{"code":"return tick.toFixed(2)+&#x27; [L/L]&#x27;"},"id":"1133","type":"FuncTickFormatter"},{"attributes":{},"id":"1061","type":"AllLabels"},{"attributes":{"bar_color":"rgba(60, 60, 60, 0.6)","end":-5.298317366548036,"format":{"id":"1137"},"js_property_callbacks":{"change:value":[{"id":"1160"}]},"sizing_mode":"stretch_width","start":-13.815510557964274,"step":0.08517193191416236,"value":-9.556913962256155},"id":"1138","type":"Slider"},{"attributes":{"code":"return tick.toFixed(2)+&#x27; [L/L]&#x27;"},"id":"1135","type":"FuncTickFormatter"},{"attributes":{"bar_color":"rgba(60, 60, 60, 0.6)","end":0.5,"format":{"id":"1139"},"js_property_callbacks":{"change:value":[{"id":"1160"}]},"sizing_mode":"stretch_width","start":-0.5,"step":0.02,"value":0},"id":"1140","type":"Slider"},{"attributes":{"code":"return Math.exp(tick).toExponential(1).toString() + &#x27; [L^2/T]&#x27;"},"id":"1137","type":"FuncTickFormatter"},{"attributes":{},"id":"1062","type":"BasicTickFormatter"},{"attributes":{"bar_color":"rgba(60, 60, 60, 0.6)","end":3.0,"format":{"id":"1141"},"js_property_callbacks":{"change:value":[{"id":"1160"}]},"sizing_mode":"stretch_width","start":0.0,"step":0.03,"value":0.0},"id":"1142","type":"Slider"},{"attributes":{"code":"return (10**tick).toFixed(2)+&#x27; [-]&#x27;"},"id":"1139","type":"FuncTickFormatter"},{"attributes":{"code":"return tick.toFixed(2)+&#x27; [-]&#x27;"},"id":"1141","type":"FuncTickFormatter"},{"attributes":{"sizing_mode":"stretch_width","style":{"font-size":"115%"},"text":"&lt;i&gt;T&lt;sub&gt;x&lt;/sub&gt;&lt;/i&gt;/&lt;i&gt;T&lt;sub&gt;y&lt;/sub&gt;&lt;/i&gt;"},"id":"1149","type":"Div"},{"attributes":{"code":"return tick.toFixed(1)+&#x27; [L]&#x27;"},"id":"1143","type":"FuncTickFormatter"},{"attributes":{"args":{"hHE":{"id":"1097"},"hHS":{"id":"1092"},"hIsoH":{"id":"1102"},"hIsoP":{"id":"1087"}},"code":"\\n    var act         = this.active\\n    hIsoH.visible   = act.includes(0)\\n    hIsoP.visible   = act.includes(1)\\n    hHE.visible     = act.includes(2)\\n    hHS.visible     = act.includes(3)\\n"},"id":"1157","type":"CustomJS"},{"attributes":{"sizing_mode":"stretch_width","style":{"font-size":"115%"},"text":"(&lt;i&gt;T&lt;sub&gt;x&lt;/sub&gt;&lt;/i&gt;\\u00b7&lt;i&gt;T&lt;sub&gt;y&lt;/sub&gt;&lt;/i&gt;)&lt;sup&gt;1/2&lt;/sup&gt;"},"id":"1148","type":"Div"},{"attributes":{"args":{"divResult":{"id":"1158"},"rG":{"id":"1153"},"sliAni":{"id":"1140"},"sliI":{"id":"1132"},"sliL":{"id":"1130"},"sliLrat":{"id":"1134"},"sliPor":{"id":"1144"},"sliQnor":{"id":"1142"},"sliT":{"id":"1138"},"sliWrat":{"id":"1136"},"srcHE":{"id":"1078"},"srcHS":{"id":"1079"},"srcIsoH":{"id":"1081"},"srcIsoP":{"id":"1082"},"srcTTD":{"id":"1080"}},"code":"// retrieve current settings\\nconst L     = sliL.value\\nconst Lrat  = sliLrat.value\\nconst Wrat  = sliWrat.value\\nconst wMax  = Lrat*L\\nconst wMin  = wMax*Wrat\\nconst Tmean = Math.exp(sliT.value)\\nconst aniso = sliAni.value\\nconst Tx    = Math.sqrt(Tmean**2 * (10**aniso))\\nconst Ty    = Math.sqrt(Tmean**2 / (10**aniso))\\nconst h2    = 0\\nconst h1    = h2 + L*Math.exp(sliI.value)\\nconst Q0    = (h1-h2)/L*Tx*(wMax-wMin)\\nconst q     = sliQnor.value*Q0/L\\nconst shape = rG.active\\nconst por   = sliPor.value\\n\\n// determine coefficients A\\nconst N = 10\\nconst M = 25\\nvar A   = getA(N,M,L,wMin,wMax,shape,Tx,Ty,h1,h2,q)\\n\\n// get fields of h and psi\\nconst [hField, pField, pMin, pMax] = getFields(L,wMin,wMax,shape,h1,h2,A,Tx,Ty);\\n\\n// extract contour line sets\\nvar lvlsH = Bokeh.LinAlg.linspace(h2,h1,11)\\nlvlsH.shift()\\nlvlsH.pop()\\n\\nconst [XsH,YsH] = getIsolines(hField,lvlsH,L,wMin,wMax,shape)\\nsrcIsoH.data[&#x27;xs&#x27;]  = XsH\\nsrcIsoH.data[&#x27;ys&#x27;]  = YsH\\nfor (var i = 0; i &lt; lvlsH.length; i++) {\\n    lvlsH[i] = (lvlsH[i]-h2)/(h1-h2) \\n} \\nsrcIsoH.data[&#x27;clr&#x27;] = lvlsH\\n\\nvar lvlsP = Bokeh.LinAlg.linspace(pMin,pMax,11)\\nlvlsP.shift()\\nlvlsP.pop()\\nconst [XsP,YsP] = getIsolines(pField,lvlsP,L,wMin,wMax,shape)\\nsrcIsoP.data[&#x27;xs&#x27;] = XsP\\nsrcIsoP.data[&#x27;ys&#x27;] = YsP\\n\\n// extract single isoline for exchange zone (he) \\nconst [heX,heY] = getIsoline(pField,pField[0][0],L,wMin,wMax,shape)\\n    \\nif (heX.length==0) {\\n    srcHE.data[&#x27;x&#x27;] = [0,L]\\n    srcHE.data[&#x27;yBot&#x27;]    = [0,0]\\n    srcHE.data[&#x27;yTop&#x27;]    = [0,0]\\n} else {\\n    srcHE.data[&#x27;x&#x27;]       = heX\\n    var heYbot = []\\n    for (var i = 0; i &lt; heX.length; i++) {\\n       heYbot.push(0) \\n    } \\n    srcHE.data[&#x27;yBot&#x27;]    = heYbot\\n    srcHE.data[&#x27;yTop&#x27;]    = heY\\n}\\n\\n// extract single isoline for hillslope zone (hs)\\nconst pVal = pField[pField.length-1][0]\\nconst [hsX,hsY] = getIsoline(pField,pVal,L,wMin,wMax,shape)\\n\\nif (hsX.length==0) {\\n    srcHS.data[&#x27;x&#x27;] = [0,L]\\n    srcHS.data[&#x27;yBot&#x27;]    = [0,0]\\n    srcHS.data[&#x27;yTop&#x27;]    = [0,0]\\n} else {\\n    const xMax = Math.max(...hsX)\\n    if ((xMax &lt; L) &amp;&amp; (Math.min(...hsY) &lt; 0.001*wMin))  {\\n        for (var i = Math.ceil(20*xMax/L); i &lt; 21; i++) {\\n            hsX.push(L*i/20)\\n            hsY.push(0)\\n        } \\n    }\\n    srcHS.data[&#x27;x&#x27;]       = hsX\\n    var hsYtop = []\\n    for (var i = 0; i &lt; hsX.length; i++) {\\n       hsYtop.push(fNorth(hsX[i],L,wMin,wMax,shape)) \\n    } \\n    srcHS.data[&#x27;yBot&#x27;]    = hsY\\n    srcHS.data[&#x27;yTop&#x27;]    = hsYtop\\n}\\n\\n// determine exchange zone area\\nvar xtemp1      = Array.from(srcHE.data[&#x27;x&#x27;])\\nvar ytemp1      = Array.from(srcHE.data[&#x27;yBot&#x27;])\\nconst xtemp2    = xtemp1.slice()\\nconst ytemp2    = Array.from(srcHE.data[&#x27;yTop&#x27;])\\nconst xPoly     = xtemp1.concat(xtemp2.reverse())\\nconst yPoly     = ytemp1.concat(ytemp2.reverse())\\nvar Aex         = polyarea(xPoly,yPoly)\\n\\n// define boolean to tell if zone exists\\nvar Qex = getQex(h1,h2,L,wMin,wMax,Tx,Ty,A)\\nvar zoneExists =  (Qex &gt;= 1e-19)\\n\\n// if zone exists: get TTD and adjust results text\\nif (zoneExists) {\\n    const [t,f] = getTTD(pField,L,wMin,wMax,h1,h2,Tx,Ty,A,shape,por);\\n\\n    // normalize ttD\\n    const Tmax = Math.max(...t)\\n    for (var i = 0; i &lt; t.length; i++) {\\n        t[i] = 100*t[i]/Tmax\\n        f[i] = 100*f[i]\\n    } \\n    srcTTD.data.t = t\\n    srcTTD.data.f = f\\n    \\n    const postQ = Math.floor(Math.log10(Qex))\\n    const preQ = Qex/(10**postQ)\\n    const postA = Math.floor(Math.log10(Aex))\\n    const preA = Aex/(10**postA)\\n    const postT = Math.floor(Math.log10(Tmax))\\n    const preT = Tmax/(10**postT)\\n    divResult.text = &#x27;The exchange flux is Q&lt;sub&gt;ex&lt;/sub&gt;  = &#x27;\\n                         + preQ.toFixed(2) \\n                         + &#x27;\\u00b710&lt;sup&gt;&#x27; \\n                         + postQ.toFixed(0) \\n                         + &#x27; &lt;/sup&gt; L&lt;sup&gt;3&lt;/sup&gt;/T. &#x27;\\n                         + &#x27;The exchange zone area is A&lt;sub&gt;ex&lt;/sub&gt; = &#x27;\\n                         + preA.toFixed(2)\\n                         + \\"\\u00b710&lt;sup&gt;\\" \\n                         + postA.toFixed(0) \\n                         +  \\"&lt;/sup&gt;  L&lt;sup&gt;2&lt;/sup&gt;.\\"\\n                         + &#x27;The longest exchange travel time is t&lt;sub&gt;max&lt;/sub&gt; = &#x27;\\n                         + preT.toFixed(2)\\n                         + \\"\\u00b710&lt;sup&gt;\\" \\n                         + postT.toFixed(0) \\n                         +  \\"&lt;/sup&gt;  T.\\"\\n} else {\\n    divResult.text = &#x27;There is no exchange zone.&#x27;\\n    srcTTD.data.t = [0,50,100]\\n    srcTTD.data.f = [0,NaN,100]\\n}\\n\\n// communicate updates\\nsrcIsoH.change.emit();\\nsrcIsoP.change.emit();\\nsrcTTD.change.emit();\\nsrcHE.change.emit();\\nsrcHS.change.emit();"},"id":"1160","type":"CustomJS"},{"attributes":{"active":[0,1,2,3],"js_property_callbacks":{"change:active":[{"id":"1157"}]},"labels":["heads","stream","zone south","zone north"],"margin":[10,3,10,3],"sizing_mode":"stretch_width"},"id":"1156","type":"CheckboxButtonGroup"},{"attributes":{},"id":"1216","type":"BasicTickFormatter"},{"attributes":{"css_classes":["divResult"],"sizing_mode":"stretch_width","style":{"font-size":"135%"},"text":"There is no exchange zone."},"id":"1158","type":"Div"},{"attributes":{},"id":"1005","type":"DataRange1d"},{"attributes":{"args":{"rG":{"id":"1153"},"sliL":{"id":"1130"},"sliLrat":{"id":"1134"},"sliWrat":{"id":"1136"},"srcLin":{"id":"1076"},"srcNorth":{"id":"1077"}},"code":"var yNorth = srcNorth.data[&#x27;y&#x27;]\\nvar xNorth = srcNorth.data[&#x27;x&#x27;]\\nconst L = sliL.value\\nconst Lrat = sliLrat.value\\nconst Wrat = sliWrat.value\\nconst wMax = Lrat*L\\nconst wMin = wMax*Wrat\\nconst val = rG.active\\n\\nfor (var i = 0; i &lt; xNorth.length; i++) {\\n    xNorth[i] = i/(xNorth.length-1)*L\\n    yNorth[i] = fNorth(xNorth[i],L,wMin,wMax,val)\\n}\\n\\nsrcLin.data[&#x27;yW&#x27;][1] = wMin\\nsrcLin.data[&#x27;yE&#x27;][1] = wMin\\nsrcLin.data[&#x27;xE&#x27;][0] = L\\nsrcLin.data[&#x27;xE&#x27;][1] = L\\nsrcLin.data[&#x27;xS&#x27;][1] = L\\n\\nsrcLin.change.emit();\\nsrcNorth.change.emit();"},"id":"1159","type":"CustomJS"},{"attributes":{},"id":"1173","type":"Selection"},{"attributes":{"children":[{"id":"1145"},{"id":"1130"}],"sizing_mode":"stretch_width"},"id":"1161","type":"Row"},{"attributes":{},"id":"1174","type":"UnionRenderers"},{"attributes":{},"id":"1175","type":"Selection"},{"attributes":{"children":[{"id":"1151"},{"id":"1134"}],"sizing_mode":"stretch_width"},"id":"1162","type":"Row"},{"attributes":{},"id":"1176","type":"UnionRenderers"},{"attributes":{},"id":"1177","type":"Selection"},{"attributes":{},"id":"1178","type":"UnionRenderers"},{"attributes":{},"id":"1179","type":"Selection"},{"attributes":{},"id":"1048","type":"WheelZoomTool"},{"attributes":{},"id":"1012","type":"BasicTicker"},{"attributes":{"axis":{"id":"1011"},"grid_line_color":null,"ticker":null},"id":"1014","type":"Grid"},{"attributes":{},"id":"1180","type":"UnionRenderers"},{"attributes":{},"id":"1022","type":"WheelZoomTool"},{"attributes":{},"id":"1181","type":"Selection"},{"attributes":{"axis_label":"y","axis_label_text_font_size":"25px","formatter":{"id":"1059"},"major_label_policy":{"id":"1058"},"major_label_text_font_size":"25px","ticker":{"id":"1016"}},"id":"1015","type":"LinearAxis"},{"attributes":{"axis":{"id":"1015"},"dimension":1,"grid_line_color":null,"ticker":null},"id":"1018","type":"Grid"},{"attributes":{"children":[{"id":"1146"},{"id":"1136"}],"sizing_mode":"stretch_width"},"id":"1163","type":"Row"},{"attributes":{},"id":"1016","type":"BasicTicker"},{"attributes":{},"id":"1182","type":"UnionRenderers"},{"attributes":{},"id":"1021","type":"SaveTool"},{"attributes":{"children":[{"id":"1147"},{"id":"1132"}],"sizing_mode":"stretch_width"},"id":"1164","type":"Row"},{"attributes":{},"id":"1183","type":"Selection"},{"attributes":{"end":100},"id":"1031","type":"Range1d"},{"attributes":{"children":[{"id":"1148"},{"id":"1138"}],"sizing_mode":"stretch_width"},"id":"1165","type":"Row"},{"attributes":{"axis":{"id":"1041"},"dimension":1,"ticker":null},"id":"1044","type":"Grid"},{"attributes":{},"id":"1019","type":"PanTool"},{"attributes":{"children":[{"id":"1149"},{"id":"1140"}],"sizing_mode":"stretch_width"},"id":"1166","type":"Row"},{"attributes":{},"id":"1184","type":"UnionRenderers"},{"attributes":{},"id":"1020","type":"ResetTool"},{"attributes":{},"id":"1185","type":"Selection"},{"attributes":{"children":[{"id":"1150"},{"id":"1142"}],"sizing_mode":"stretch_width"},"id":"1167","type":"Row"},{"attributes":{},"id":"1047","type":"SaveTool"},{"attributes":{},"id":"1042","type":"BasicTicker"},{"attributes":{"children":[{"id":"1152"},{"id":"1144"}],"sizing_mode":"stretch_width"},"id":"1168","type":"Row"},{"attributes":{"active_drag":null,"active_multi":null,"active_scroll":{"id":"1022"},"logo":"grey","tools":[{"id":"1045"},{"id":"1046"},{"id":"1047"},{"id":"1048"}]},"id":"1049","type":"Toolbar"},{"attributes":{},"id":"1186","type":"UnionRenderers"},{"attributes":{},"id":"1203","type":"Title"},{"attributes":{"children":[{"id":"1161"},{"id":"1162"},{"id":"1163"},{"id":"1164"},{"id":"1165"},{"id":"1166"},{"id":"1167"},{"id":"1168"},{"id":"1153"},{"id":"1154"},{"id":"1156"},{"id":"1158"}],"sizing_mode":"stretch_width"},"id":"1169","type":"Column"},{"attributes":{},"id":"1045","type":"PanTool"},{"attributes":{"fill_alpha":0.1,"fill_color":"rgb(201, 174, 105)","x":{"field":"x"},"y1":{"field":"yTop"},"y2":{"field":"yBot"}},"id":"1091","type":"VArea"},{"attributes":{"tabs":[{"id":"1170"},{"id":"1171"}]},"id":"1172","type":"Tabs"},{"attributes":{"child":{"id":"1002"},"title":"Flow Net"},"id":"1170","type":"Panel"},{"attributes":{"data_source":{"id":"1082"},"glyph":{"id":"1085"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"1086"},"view":{"id":"1088"}},"id":"1087","type":"GlyphRenderer"},{"attributes":{"child":{"id":"1028"},"title":"Travel Times"},"id":"1171","type":"Panel"},{"attributes":{"axis_label":"F(t) in %","axis_label_text_font_size":"25px","formatter":{"id":"1213"},"major_label_policy":{"id":"1212"},"major_label_text_font_size":"25px","ticker":{"id":"1042"}},"id":"1041","type":"LinearAxis"},{"attributes":{},"id":"1046","type":"ResetTool"},{"attributes":{"source":{"id":"1082"}},"id":"1088","type":"CDSView"},{"attributes":{},"id":"1007","type":"LinearScale"},{"attributes":{"active_drag":null,"active_multi":null,"active_scroll":{"id":"1022"},"logo":"grey","tools":[{"id":"1019"},{"id":"1020"},{"id":"1021"},{"id":"1022"}]},"id":"1023","type":"Toolbar"},{"attributes":{},"id":"1212","type":"AllLabels"},{"attributes":{"sizing_mode":"stretch_width","style":{"font-size":"115%"},"text":"&lt;i&gt;w&lt;sub&gt;max&lt;/sub&gt;&lt;/i&gt;/&lt;i&gt;w&lt;sub&gt;min&lt;/sub&gt;&lt;/i&gt;"},"id":"1146","type":"Div"},{"attributes":{},"id":"1213","type":"BasicTickFormatter"},{"attributes":{},"id":"1038","type":"BasicTicker"},{"attributes":{"end":100},"id":"1029","type":"Range1d"},{"attributes":{},"id":"1033","type":"LinearScale"},{"attributes":{},"id":"1035","type":"LinearScale"},{"attributes":{"axis_label":"t in % of tMax","axis_label_text_font_size":"25px","formatter":{"id":"1216"},"major_label_policy":{"id":"1215"},"major_label_text_font_size":"25px","ticker":{"id":"1038"}},"id":"1037","type":"LinearAxis"},{"attributes":{"line_alpha":{"value":0.1},"line_color":"rgba(230,232,237, 1.0)","line_width":{"value":5},"xs":{"field":"xs"},"ys":{"field":"ys"}},"id":"1086","type":"MultiLine"},{"attributes":{"data_source":{"id":"1079"},"glyph":{"id":"1090"},"hover_glyph":null,"muted_glyph":null,"nonselection_glyph":{"id":"1091"},"view":{"id":"1093"}},"id":"1092","type":"GlyphRenderer"},{"attributes":{"axis":{"id":"1037"},"ticker":null},"id":"1040","type":"Grid"},{"attributes":{"sizing_mode":"stretch_width","style":{"font-size":"115%"},"text":"&lt;i&gt;\\u0394h/L&lt;/i&gt;"},"id":"1147","type":"Div"},{"attributes":{},"id":"1215","type":"AllLabels"},{"attributes":{"fill_alpha":0.4,"fill_color":"rgb(201, 174, 105)","x":{"field":"x"},"y1":{"field":"yTop"},"y2":{"field":"yBot"}},"id":"1090","type":"VArea"},{"attributes":{"source":{"id":"1079"}},"id":"1093","type":"CDSView"},{"attributes":{"line_alpha":{"value":0.7},"line_color":"rgba(230,232,237, 1.0)","line_width":{"value":5},"xs":{"field":"xs"},"ys":{"field":"ys"}},"id":"1085","type":"MultiLine"}],"root_ids":["1169","1172"]},"title":"Bokeh Application","version":"2.3.2"}}';
                var render_items = [{"docid":"aad0f763-9f6d-4b1d-a19f-82e554a0c78f","root_ids":["1169","1172"],"roots":{"1169":"9f8985b9-7bb5-489d-b4e8-55abac1ebc87","1172":"56f15a42-282b-43f8-9854-d15ae355da39"}}];
                root.Bokeh.embed.embed_items(docs_json, render_items);
              
                }
                if (root.Bokeh !== undefined) {
                  embed_document(root);
                } else {
                  var attempts = 0;
                  var timer = setInterval(function(root) {
                    if (root.Bokeh !== undefined) {
                      clearInterval(timer);
                      embed_document(root);
                    } else {
                      attempts++;
                      if (attempts > 100) {
                        clearInterval(timer);
                        console.log("Bokeh: ERROR: Unable to run BokehJS code because BokehJS library is missing");
                      }
                    }
                  }, 10, root)
                }
              })(window);
            });
          };
          if (document.readyState != "loading") fn();
          else document.addEventListener("DOMContentLoaded", fn);
        })();