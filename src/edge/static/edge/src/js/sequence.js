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
              return (row_start >= a.base_first && row_start <= a.base_last) ||
                     (row_end >= a.base_first && row_end <= a.base_last) ||
                     (a.base_first >= row_start && a.base_first <= row_end);
            });

            _.map(row_annotations, function(a) {
              var tr = $('<tr></tr>');
              var aspan = row.length;
              if (a.base_first > row_start) {
                var td = $('<td colspan="'+(a.base_first-row_start)+'"></td>');
                tr.append(td);
                aspan -= (a.base_first-row_start);
              }
              if (a.base_last < row_end) {
                aspan -= (row_end-a.base_last);
              }
              var td = $('<td colspan="'+aspan+'" class="sequence-annotation"></td>');
              td.append(a.display_name);
              tr.append(td);
              if (a.base_last < row_end) {
                var td = $('<td colspan="'+(row_end-a.base_last)+'"></td>');
                tr.append(td);
              }
              $('table', selector).append(tr);
            });

            var tr = $('<tr></tr>');
            _.map(row.split(''), function(bp) {
                var td = $('<td class="sequence-data"></td>');
                td.append(bp);
                tr.append(td);
            });
            $('table', selector).append(tr);

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
