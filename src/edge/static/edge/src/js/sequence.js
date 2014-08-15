window.SequenceViewer = function ($, options) {
    var dom_id = 'sequence-viewer';
    if ('dom_id' in options) { dom_id = options.dom_id; }
    var selector = '#'+dom_id;
    $(selector).html('<table></table>');

    function clear() {
        $('table', selector).empty();
    }

    function setSequenceWithAnnotations(sequence, annotations, start_bp, rowlen) {
        rowlen = typeof rowlen !== 'undefined' ? rowlen : 80;
        clear();
        var rsplit = new RegExp('(.{1,'+rowlen+'})', 'g');
        var rows = sequence.match(rsplit);

        var row_start = start_bp;
        _.map(rows, function(row) {
            var row_end = row_start+row.length-1;
            var row_annotations = _.select(annotations, function(a) {
                return (row_start <= a.base_first && a.base_first <= row_end) ||
                       (row_start <= a.base_last && a.base_last <= row_end);
            });

            var start_annotations = [];
            var end_annotations = [];

            _.map(row_annotations, function(a) {
                if (row_start <= a.base_first && a.base_first <= row_end) {
                    var i = a.base_first-row_start;
                    if (start_annotations[i] === undefined) { start_annotations[i] = []; }
                    start_annotations[i].push(a);
                }
                if (row_start <= a.base_last && a.base_last <= row_end) {
                    var i = a.base_last-row_start;
                    if (end_annotations[i] === undefined) { end_annotations[i] = []; }
                    end_annotations[i].push(a);
                }
            });

            var tr_top = $('<tr class="sequence-annotation sequence-annotation-start"></tr>');
            var tr_mid = $('<tr class="sequence-bp"></tr>');
            var tr_bot = $('<tr class="sequence-annotation sequence-annotation-end"></tr>');

            var i = 0;
            _.map(row.split(''), function(bp) {
                var td = $('<td></td>');
                td.append(bp);
                tr_mid.append(td);

                if (start_annotations.length > 0) {
                  var td = $('<td></td>');
                  if (start_annotations[i]) {
                    _.map(start_annotations[i], function(a) {
                      var a = $('<div></div>').append('['+a.display_name);
                      td.append(a);
                    });
                  }
                  tr_top.append(td);
                }

                if (end_annotations.length > 0) {
                  var td = $('<td></td>');
                  if (end_annotations[i]) {
                    _.map(end_annotations[i], function(a) {
                      var a = $('<div></div>').append(']'+a.display_name);
                      td.append(a);
                    });
                  }
                  tr_bot.append(td);
                }

                i += 1;
            });

            if (start_annotations.length > 0) { $('table', selector).append(tr_top); }
            $('table', selector).append(tr_mid);
            if (end_annotations.length > 0) { $('table', selector).append(tr_bot); }

            row_start += row.length;
        });
    }

    function setSequence(sequence, start_bp, rowlen, column) {
        rowlen = typeof rowlen !== 'undefined' ? rowlen : 80;
        column = typeof column !== 'undefined' ? column : 10;
        clear();

        var rsplit = new RegExp('(.{1,'+rowlen+'})', 'g');
        var csplit = new RegExp('(.{1,'+column+'})', 'g');
        var rows = sequence.match(rsplit);
        var tr = $('<tr></tr>');

        var s = start_bp;
        var td_s = $('<td class="sequence-pos"></td>');
        var td_d = $('<td class="sequence-data"></td>');
        var td_l = $('<td class="sequence-pos"></td>');
        for (var i=0; i<rows.length; i++) {
            var l = s+rows[i].length-1;
            td_s.append(s+'/'+(s-start_bp+1)+'<br/>');
            var columns = rows[i].match(csplit);
            for (var j=0; j<columns.length; j++) {
                td_d.append('<span>'+columns[j]+'</span>');
            }
            td_d.append('<br/>');
            td_l.append(l+'/'+(l-start_bp+1)+'<br/>');
            s = l+1;
        }
        tr.append(td_s);
        tr.append(td_d);
        tr.append(td_l);
        $('table', selector).append(tr);
    }

    return {
        clear: clear,
        setSequence: setSequence,
        setSequenceWithAnnotations: setSequenceWithAnnotations
    }
}
