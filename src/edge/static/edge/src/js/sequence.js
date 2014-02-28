window.SequenceViewer = function ($, options) {
    var dom_id = 'sequence-viewer';
    if ('dom_id' in options) { dom_id = options.dom_id; }
    var selector = '#'+dom_id;
    $(selector).html('<table></table>');

    function clear() {
        $('table', selector).empty();
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
        setSequence: setSequence
    }
}
