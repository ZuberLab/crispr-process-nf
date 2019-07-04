#!/usr/bin/env python3

"""fileutilities.py

Author: Kimon Froussios
Compatibility tested: python 3.5.2
Last reviewed: 10/05/2019

This module is a solution for Frequently Performed Generic Tasks that involve
multiple files:
* repeating a command for a range of files (not currently parallelized),
* accessing and restructuring (multiple) delimited files.
* miscellaneous stuff. Some of it is auxiliary to the primary functions, some
  is a legacy of this module's evolution of concept.

The module provides a library of flexible functions as
well as a main() implementing the primary use scenarios.

Execute with -h in a shell to obtain syntax and help.
"""

# This module consists of:
# -    a class for handling lists of files,
# -    a library of functions that perform generic tasks on multiple files, and
# -    a main that provides access to most of the above functionality

# NOTE about DataFrame indexes and headers:
# Although dataframes support row and column labels, these make content manipulations
# in the context of this module harder. Instead, any labels present in the text
# input are treated as plain rows or columns. These may be optionally dropped or
# preserved, but the dataframe in-built labels for columns and rows are reserved
# solely for custom use. When appropriate these custom labels will be included in
# the output.


import os, sys, string, re, subprocess, random, argparse
import pandas as pd
from builtins import list
from collections import Counter

import mylogs as ml


#####   F U N C T I O N S   #####


# http://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
def natural_sorted(l):
    """Sort list of numbers/strings in human-friendly order.

    Args:
        l(list): A list of strings.
    Returns:
        list
    """
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


def expand_fpaths(flist):
        """Fully expand and absolute-ify the paths of listed files.

        Does not verify path validity. All paths are expanded.

        Args:
            flist[str]: A list/FilesList of files.
        Returns:
            [str]: List of expanded paths.
        """
        return [os.path.abspath(os.path.expanduser(str(f))) for f in flist]


# Helper function, string check.
def istext(s):
    """Check if a string is (probably) text.

    Use heuristic based on characters contained in the string, adapted from:
    http://code.activestate.com/recipes/173220-test-if-a-file-or-string-is-text-or-binary/

    Args:
       s(str): A string to test.
    Returns:
        bool
    """
    # Copy-pasted. No idea what the code means.
    text_characters = "".join(list(map(chr, list(range(32, 127)))) + list("\n\r\t\b"))
    _null_trans = string.maketrans("", "")
    if "\0" in s:
       return False
    if not s:  # Empty files/strings are considered text
        return True
    # Get the non-text characters (maps a character to itself then
    # use the 'remove' option to get rid of the text characters.)
    t = s.translate(_null_trans, text_characters)
    # If more than 30% non-text characters, then
    # this is considered a binary file
    if float(len(t))/float(len(s)) > 0.30:
        return False
    return True


def slink(flist, aliases=None, dir="./", autoext=True):
    """Create symbolic links for multiple files.

    Create a link for each of the listed paths into the specified directory,
    using the specified aliases. Items in the lists will be matched one for
    one.
    If the aliases argument is omitted, the names for the links will be drawn
    from the aliases attribute of the paths list, if it is a FilesList object.
    If no aliases exist in either form, the files will be linked in the current
    or specified directory, using names their current basename.

    If linking to files of the same name located in different directories, a
    number will be automatically suffixed to the basename.

    Args:
        flist[str]: A list/FilesList of paths to link to.
        aliases[str]: A list of respective names for the created links. If
                    omitted, the alias attribute of the flist argument will be
                    used, and failing that, the existing basenames will be used.
        dir(str): The path to the directory in which the links should be
                    placed. (Default "./")
        autoext(bool): Add the file extensions to the created links, if the
                    links are created from aliases that lack them.
                    (Default True)
    """
    if not aliases:
        # No link names provided. Try to find them elsewhere or create them.
        try:
            # flist is a FilesList and has the aliases attribute.
            aliases = flist.aliases
        except AttributeError:
            # flist is a plain list, so compute the link name from the file name.
            aliases = [os.path.basename(p) for p in flist]
    # Check for duplicate aliases and amend them.
    # This applies mainly to link names automatically created from filenames, as
    # the same file name can exist in different directories.
    if len(set(aliases)) < len(flist):
        aliases = autonumerate(aliases)
    # Add extensions where necessary, if desired.
    if autoext:
        for i in range(0, len(flist)):
            (b, p) = os.path.splitext(flist[i])
            c = p
            # If it's a .gz, include the next nested extension as well.
            if p == ".gz":
                p = os.path.splitext(b)[1] + p
            # Don't duplicate the extension if the alias already has it.
            a = os.path.splitext(aliases[i])[1]
            if c != a:
                aliases[i] = aliases[i] + p
    # Link.
    for i, mypath in enumerate(flist):
        os.symlink(mypath, os.path.join(dir, aliases[i]))


# Helper function.
def autonumerate(things):
    """Detect duplicate entries in a string list and suffix them.

    Suffixes are in _N format where N a natural number >=2. Existing suffixes
    in that format will also be detected and incremented.

    Args:
        things[str]: A list of strings.
    Returns:
        [str]: A corrected list of strings.
    """
    c = Counter(things);
    # Because I use decrement, reversing the list ensures first instance gets smallest number.
    things.reverse()
    for i, t in enumerate(things):
        n = c[t]
        if n > 1:  # The first occurrence is not suffixed.
            newname = t +'_' + str(n)
            while newname in things:  # Check for already present suffixes
                n += 1
                newname = t +'_' + str(n)
            things[i] = newname
            c[t] -= 1
    things.reverse()
    return things


def make_names(items, parameters):
    """Automatically create file names based on parameters.

    If automatic names happen to turn out identical with one another, unique
    numbers are appended to differentiate them. Check documentation for
    autonumerate().

    Args:
        items[str]: A list of strings/filenames/paths to use as the basis for
                    the output names.
        parameters(str,str,str): The first element is the output directory,
                    the second is a common prefix to add to the names,
                    the third is a common suffix to add to the names.
                    Like so: <out[0]>/<out[1]>item<out[2] .
                    If any of the 3 values in None, no outnames will be made.
                    Use current directory and empty strings as necessary.
    Returns:
        [str]: A list of file paths.
    """
    outfiles = []
    if None not in parameters:
        for i in items:
            outfiles.append(os.path.join(os.path.abspath(os.path.expanduser(parameters[0])),
                                         parameters[1] + i + parameters[2]) )
        autonumerate(outfiles)
    return outfiles


def do_foreach(flist, comm, comments=False, progress=True, out=(None,None,None), log=False):
    """Execute an arbitrary command for each of the listed files.

    Enables executing a shell command over a range of items, by inserting the
    item values into the command as directed by place-holder substrings.
    Although the above is how it is meant to be used, the values in the
    FilesList could be overridden to be any arbitrary string, in which case,
    only {val} will have the desired effect. The other placeholder values are
    computed with the assumption of the values being files, so may not be
    sensible when the items are not files.

    This is the only function with comments or progress attributes, because
    invoked commands print their own output directly, so any informative messages
    controlled by this library will need to be inserted in real time.

    Args:
        flist[]: A FilesList.
        comm[str]: The components of an arbitrary command, with place-holders:
                    {abs} : absolute path of file.
                    {dir} : absolute path of the file's directory.
                    {val} : the actual value specified as target
                    {bas} : the basename of the file, without the last extension.
                    {alias}: the alias for the file, if iterating through a FilesList.
                    Placeholders can be nested, to allow nested calls of fileutilities:
                    i.e. {{abs}}. A layer of nesting is peeled off each time the function is called,
                    until the placeholders are allowed to be evaluated.
        comments(bool): Print commented call details to STDOUT. (Default False)
        progress(bool): Show start and completion of iterations on STDERR.
                    (Default True)
        out(str,str,str): The first element is the output directory, the second
                    is a common prefix to add to the names, the third is a
                    common suffix to add to the names. Check documentation for
                    make_names().
        log(bool): Log to /commands.log each individual call.
    """
    outstream= sys.stdout
    # Create output files. [] if out contains None.
    outfiles = make_names(flist, out)
    for i, (myfile, myalias) in flist.enum():
        # Substitute place-holders.
        command = []
        for c in comm:
            # Evaluate placeholders, if they are not nested.
            (mypath, mybase) = os.path.split(str(myfile))
            c = re.sub(r"(?<!\{){abs}(?!\})", str(myfile), c)
            c = re.sub(r"(?<!\{){dir}(?!\})", mypath, c)
            c = re.sub(r"(?<!\{){val}(?!\})", mybase, c)
            c = re.sub(r"(?<!\{){bas}(?!\})", os.path.splitext(mybase)[0], c)
            c = re.sub(r"(?<!\{){ali}(?!\})", str(myalias), c)
            # Peel off a layer of nesting for the remaining placeholders and flags.
            c = c.replace('{{abs}}', '{abs}')
            c = c.replace('{{dir}}', '{dir}')
            c = c.replace('{{val}}', '{val}')
            c = c.replace('{{bas}}', '{bas}')
            c = c.replace('{{ali}}', '{ali}')
            c = c.replace(',-', '-')
            # This argument is ready to go now.
            command.append(c)
        # Redirect output.
        if outfiles:
            outstream = open(outfiles[i], 'w')
        # Verbose stuff.
        try:
            see = " ".join(command)
            if log:
                ml.log_message(message=see, logfile="./subcommands.log")
            if comments and out == (None,None):
                outstream.write(ml.infostring("CWD: "+ os.getcwd() +"\tDO: "+ see))
            if progress:
                sys.stderr.write(ml.infostring("DO: "+ see))
        except IOError:
            pass
        # Do the thing.
        subprocess.call(" ".join(command), stdout=outstream, shell=True)
        # Optionally identify iteration.
        try:
            if comments and out == (None,None):
                outstream.write(ml.infostring("Finished: "+ str(myalias) +"\n"))
            if progress:
                sys.stderr.write(ml.infostring("Finished: "+ str(myalias) +"\n"))
        except IOError:
            pass
        finally:
            if outfiles:
                outstream.close()


def swap_strFiles(flist, insep=[","], outsep="\t"):
    """Replace the column separator with a different one.

    Supports multiple different delimiters in the input, to support one-step
    uniformity when the input files have different delimiters, but ALL input
    will be split at ALL/ANY occurring delimiters. If the delimiter of one
    file is present in a different use in an other file, the output may not
    be what you want.
    Although made for converting delimited text, inseps and outsep could be any
    substring in a text, delimited or not.

    Args:
        flist: A list/FilesList of delimited text files.
        insep[str]: A list of regex strings. (Default [","])
        outsep(str): New column separator. (Default "\t")
    Returns:
        [str]: A list of strings with the changed delimiters. One string per
                file. It is up to you to decide  what to do with the new
                strings. The order of strings is the same as the input.
    """
    input = []
    if flist == []:
        # Read all data from STDIN at once. Input[] gets a single entry.
        input.append(sys.stdin.read())
    else:
        # Read in all the files entirely. Input[] gets as many entries as there are files.
        for myfile in flist:
            with open(myfile) as f:
                input.append(f.read())
    return swap_substr(input, insep, outsep)

# Helper function
def swap_substr(slist, insep=[","], outsep="\t"):
    """Replace all occurrences of insep with outsep.

    Insep may be a regex.

    Args:
        slist[str]: A list of strings.
        insep[str]: A list of regex strings. (Default [","])
        outsep(str): New substring. (Default "\t")
    Returns:
        [str]: A list of the edited strings. The order of the strings is the
                    same as the input.
    """
    rx = re.compile("|".join(insep), re.MULTILINE)
    result = []
    for s in slist:
        # Replace all delimiters with the new one.
        result.append(rx.sub(outsep, s))
    return result


def prepare_df(df, myalias="", keyCol=None, keyhead="row_ID", header=False, cols=None, appendNum=True):
    """Prepare row names and column names.

    Assign column as row labels, rename columns based on their position and an
    arbitrary alias name for the dataframe, drop the first row.

    Args:
        df(pandas.DataFrame): A dataframe.
        myalias(str): The basename for the relabelling.
        header(bool): Remove first row (Default False).
        keyCol(int): Column to be used as row index. If None, no index will be
                    used. (Default None)
        keyhead(str): Label for the index.
        cols[int]: Custom index numbers for the columns (Default None). If None
                    then their current index positions are used.
        appendNum(bool): Append the columns' positional indices to the alias
                    when making the new names (True).
    Returns:
        pandas.DataFrame
    """
    # Set row labels.
    if keyhead is None:
        keyhead = "row_ID"
    if keyCol is not None:
        # Add index without dropping it, so as not to affect column positions.
        df.set_index(df.columns.values.tolist()[keyCol], inplace=True, drop=False)
        df.index.name = str(keyhead)
    # Make custom column labels, based on alias and column position.
    if not cols:
        cols = list(range(0, df.shape[1]))
    labels = []
    if appendNum:
        labels = [str(myalias) +"_|"+ str(i) for i in cols]
    else:
        labels = [str(myalias) for i in cols]
    df.columns = labels
    # Remove header.
    if header:
        df.drop(df.index.values.tolist()[0], axis=0, inplace=True)
    return df


def count_columns(flist=[None], colSep=["\t"]):
    """Determine the number of fields in each file by inspecting the first row.

    Args:
        flist: A list of FilesList of files.
        colSep[str]: A list of characters used to separate columns.
    Returns:
        [int]: A list, in the same order as the given files.
    """
    tokenizer = re.compile("|".join(colSep))
    counts = []
    for file in flist:
        f = None
        if file is None:
            f = sys.stdin
            file = "<STDIN>"
        else:
            f = open(file)
        while True:
            line = f.readline()
            # Skip comments.
            if line[0] != "#":
                counts.append(len( tokenizer.split(line.rstrip()) ))
                break
            f.readline
        if f != sys.stdin:
            f.close()
    return counts


def get_valuesSet(flist=[None], axis='r', index=0, filter='a', colSep=["\t"]):
    """"List the set of different values in the column(s)/row(s).

    Args:
        flist: A list of FilesList of files.
        colSep[str]: A list of characters used to separate columns.
        index: Position index of the required column/row.
        axis(str): Data slice orientation - 'r' for row, 'c' for column.
        filter(str): non redundant set of: 'a' - all, 'u' - unique, 'r' -
                    repeated values.
    Returns:
        [[]]: A list of lists. The inner lists represent the sets in order as
                    requested.
    Raises:
        ValueError: Invalid axis or filter values.

    """
    tokenizer = "|".join(colSep)
    result = []
    if flist == []:
        # Use None as a flag to read from STDIN
        flist.append(None)
    # Test if it is a FilesList or plain list. Upgrade it if it's plain.
    # It will have at least one entry for sure by now, either way.
    try:
        flist.aliases[0]
    except AttributeError:
        flist = FilesList(flist)
    # Main part of this function.
    results = []
    for f, (myfile, myalias) in flist.enum():
        # Input.
        df = None
        instream = sys.stdin
        if myfile is not None:
            instream = open(myfile)
        df = pd.read_csv(instream, sep=tokenizer, header=None, index_col=None, comment="#", engine='python')
        if instream != sys.stdin:
            instream.close()
        # Get value set.
        values = None
        if axis == 'r':
            values = df.iloc[int(index),:].tolist()
        elif axis == 'c':
            values = df.iloc[:,int(index)].tolist()
        else:
            raise ValueError("".join(["Unrecognized option: axis=", axis]))
        # Get count per value
        c = Counter(values);
        # Filter.
        if filter == 'a':
            results.append( set(values) )
        elif filter == 'u':
            results.append( set([v for v in values if c[v] == 1]) ) # set() is redundant but keeps output type consistent
        elif filter == 'r':
            results.append( set([v for v in values if c[v] > 1]) )
        else:
            raise ValueError("".join(["Unrecognized option: filter=", filter]))
    return results


def get_columns(flist=[None], cols=[0], colSep=["\t"], header=False, index=None, merge=True):
    """Obtain the specified columns.

    Comment lines starting with '#' are ignored.
    The data columns are assembled into a single DataFrame.

    The returned columns will be labeled based on the name of the file they
    came from and their position in it. Existing labels are optionally
    preserved as the top row or can be skipped entirely.

    If an index is specified, it will be used only for merging, and will NOT be
    included in the output columns, unless explicitly present in cols[].

    Args:
        flist: A list/FilesList of delimited plain text files.
        header(bool): Crop the header line (first non-comment line). (Default False)
        cols[int/str] : A list of positional indexes or names or ranges of the
                    desired columns. (Default [0]).
        colSep[str]: List of characters used as field separators.
                    (Default ["\t"]).
        merge(bool): Concatenate results from all files into a single
                    dataframe. If False, a list of dataframes is returned
                    instead. (Default True).
        index(int): Column to be used as row index for merging. (Default None)
    Returns:
        [pandas.DataFrame]: List of DataFrames. If merge=True, only the
                    first element will be populated.
    """
    tokenizer = "|".join(colSep)
    result = []
    if flist == []:
        # Use None as a flag to read from STDIN
        flist.append(None)
    # Test if it is a FilesList or plain list. Upgrade it if it's plain.
    # It will have at least one entry for sure by now, either way.
    try:
        flist.aliases[0]
    except AttributeError:
        flist = FilesList(flist)
    # Parse.
    keyhead = None
    for f, (myfile, myalias) in flist.enum():
        # I used to use the pandas parser, with my own parser used only as fallback
        # for problematic cases. As flexibility requirements increased, using the
        # pandas parser became too opaque and difficult to maintain,
        # so now all cases are delegated to mine.
        df = get_columns_manual(myfile, cols=cols, colSep=colSep, header=header,
                                        alias=myalias, index=index)
        if not keyhead:
            keyhead = df.index.name
        result.append(df)
    # Merge.
    if merge:
        result = [pd.concat(result, axis=1, join='outer', ignore_index=False, sort=False), ]
        result[0].index.name = keyhead
    return result


# Helper function
def get_columns_manual(file=None, cols=[0], colSep=["\t"], header=False, index=None, alias=None):
    """Get specified columns from a file where rows have varying numbers of fields.

    Some tables contain a fixed set of fields followed by optional fields. In
    these rare cases, traditional parsers fail due to inconsistent number of
    fields. This function provides a work-around for that.

    It is entirely the user's responsibility to ensure that the inconsistent
    row lengths are not a symptom of table corruption/malformation and that it
    is safe and reliable to extract the desired columns. If a row is shorter
    than expected, it is padded with the value "IDXERROR". If this value shows
    up in your result and you are not explicitly expecting it, you should stop
    and seriously examine your input table.

    Args:
        file(str): A delimited plain text file.
        header(bool): If True, the first non-comment line will not be in
                    the data. (Default False)
        cols[int]: A list of positional indexes of the desired columns.
                    (Default [0]).
        colSep[str]: List of regex strings for field separators.
                    (Default ["\t"]).
        index(int): Position of column to be used as row index. (Default None)
        alias(str): An alias for the file. Used for naming the columns.
    Returns:
        pandas.DataFrame: DataFrame with the columns, labeled by original
                    column number, ordered as specified.
    """
    tokenizer = re.compile("|".join(colSep))
    # Input source.
    f = None
    if file is None:
        f = sys.stdin
        file = "STDIN"
    else:
        f = open(file)
    if alias is None:
        alias = FilesList.autoalias(file)
    # Import data.
    keyhead = None
    values = []
    labels = []
    for l, line in enumerate(f):
        if line[0] == '#' or line == "\n":
            # Skip comments and empty lines.
            continue
        else:
            # Get the fields.
            fields = tokenizer.split(line.rstrip("\n"))
            # Column labels from the first non-comment non-empty row,
            # regardless of whether they really are labels or not.
            if not labels:
                labels = fields
            # Find out name of row index.
            if (not keyhead) and header and (index is not None):
                keyhead = str(fields[index])
            # Get columns.
            selection = []
            expandedcols = []
            for c in cols:
                v = str(c).split(":")
                if len(v) == 1:
                    try:
                        expandedcols.append(int(v[0]))
                    except ValueError:
                        expandedcols.append(labels.index(v[0]))
                else:
                    try:
                        expandedcols.extend(list(range(int(v[0]), int(v[1]) + 1)))
                    except TypeError:
                        expandedcols.extend(list(range(labels.index(v[0]), labels.index(v[1]) + 1)))
            for i in expandedcols:
                try:
                    selection.append(fields[i])
                except IndexError:
                    # Silently adding fields is too dangerous, so a flag value is needed.
                    # Values like None or NA can sometimes be legitimate values for fields.
                    selection.append("IDXERROR")
            # Add the key at the end, where they won't interfere with column numbers.
            if index is not None:
                selection.append(fields[index])
            values.append(selection)
    if f != sys.stdin:
        f.close()
    # Adjust index of row keys to reflect the fact I stuck them at the end.
    if index is not None:
        index = len(values[0])-1
        expandedcols.append("my_garbage_label_row_key")
    # Package data nicely.
    df = pd.DataFrame(data=values)
    df.astype(str, copy=False)   		# Uniform string type is simplest and safest.
    df = prepare_df(df, myalias=alias, keyCol=index, header=header, cols=expandedcols,
                    keyhead=keyhead, appendNum=True if len(expandedcols)>1 else False)
    if index is not None:
        df.drop(alias+"_|my_garbage_label_row_key", 1, inplace=True)
    return df


def get_random_columns(flist, colSep=["\t"], k=1, header=False, index=None, merge=True):
    """ Get k random columns from each file.

    The returned columns will be labeled based on the name of the file they
    came from and their position in it. Existing labels are optionally
    preserved as the top row or can be skipped entirely.

    If an index is specified, it will be used for merging (if applicable) and
    will be included as a column in each output file.

    Args:
        flist: A list or FilesList of files.
        k(int): How many columns to get.
        colSep[str]: A list of characters used as field separators.
                    (Default ["\t"])
        header(bool): Strip column headers. (Default False)
        index(int): Column to use as row index for merging. (Default None)
        merge(bool): Concatenate results from all files into a single
                    dataframe. If False, a list of dataframes is returned
                    instead. (Default True).
    Returns:
        [pandas.DataFrame]: List of DataFrames. If merge=True, only the
                        first element will be populated.
    """
    tokenizer = "|".join(colSep)
    # The files may have different number of columns
    fieldNums = count_columns(flist, colSep)
    result = []
    if flist == []:
        # Use None as a flag to read from STDIN
        flist.append(None)
    keyhead = None
    # Test if it is a FilesList or plain list. Upgrade it if it's plain.
    # get_columns() does this too, but as I call it per item in flist, I *must*
    # preserve any alias that is potentially already specified.
    try:
        flist.aliases[0]
    except AttributeError:
        flist = FilesList(flist)
    # Get columns.
    for f, (myfile, myalias) in flist.enum():
        cols = []
        if index is not None:
            # Generate random choice of columns.
            # range() is right-open.
            cols.extend(random.sample(list(range(0, fieldNums[f]-1)), k))
        else:
            cols = random.sample(list(range(0,fieldNums[f])), k)
        # Would normally delegate the actual getting to get_columns() but there
        # are too many little differences to accommodate that complicate the code
        # to the point of defeating any benefits from code re-use.
        df = pd.read_csv(myfile, sep=tokenizer, header=None, index_col=None, comment="#", engine='python')
        if (not keyhead) and header and (index is not None):
            keyhead = str(df.iloc[0,index])
        # Adjust row and column labels.
        df = prepare_df(df, myalias=myalias, keyCol=index, header=header, keyhead=keyhead,
                        appendNum=True if k>1 else False)
        # Slice the part I need.
        df = df.iloc[:,cols]
        result.append(df)
    # Merge.
    if merge:
        result = [pd.concat(result, axis=1, join='outer', ignore_index=False)]
        result[0].index.name = keyhead
    return result


def append_columns(flist, colSep=["\t"], header=False, index=None, merge=True, type='outer'):
    """Append all columns from the files, as they are.

    Inner or outer concatenation by a unique index, or index-less concatenation.
    Ideally used with a unique index and for only for same-length files.

    Args:
        flist: A list/FilesList of files to combine.
        colSep[str]: A list of characters used as field delimiters.
                    (Default ["\t"])
        header(bool): First non-comment line as column labels. (Default False)
        index(int): Column to use as row index (same in all files).
                    (Default None)
                    If None, the number of rows can differ between files and will be
                    padded (outer) or truncated (inner), otherwise the row number must
                    be the same in all files.
        type(str): Join type 'inner' or 'outer'.
    Returns:
        pandas.Dataframe
    """
    try:
        flist.aliases[0]
    except AttributeError:
        flist = FilesList(flist)
    # Determine how many columns each file has.
    numofcols = count_columns(flist, colSep=colSep)
    # Delegate fetching all the columns.
    data = []
    keyhead = None
    for f, (myfile, myalias) in flist.enum():
        # List the columns and remove the index one from among them.
        cols = [i for i in range(0,numofcols[f]) if i != index]
        df =get_columns(FilesList(files=[myfile], aliases=[myalias]), cols=cols,
                     colSep=colSep, header=header, merge=False, index=index)[0]
        data.append( df )
    # Merge. Row indexes will have been assigned by get_columns(), if applicable.
    keyhead = data[0].index.name
    result = pd.concat(data, axis=1, join=type, ignore_index=False, sort=False)
    result.index.name = keyhead
    return result


def merge_tables(flist, colSep=["\t"], header=False, index=0, merge=True, type='outer', saveHeader=False, dedup=True):
    """Incrementally merge tables.

	Join the first two files and then join the third file to the merged first two, etc.
	For assymetric joins (left or right) the order of files in flist can change the outcome.
	For symmetric joins (inner or outter) the order of files should not change the outcome.

    Args:
        flist: A list/FilesList of files to combine.
        colSep[str]: A list of characters used as field delimiters.
                     (Default ["\t"])
        header(bool): Crop first non-comment line as column labels. (Default False)
        index(int): Column to use as row index (same in all files).
                    (Default 0)
        type(str): 'left', 'right', 'outer' or 'inner' merge. (Default outer)
        saveHeader(bool): Exclude the first row from sorting upon merging. (False)
                          Necessary when the header is not to be cropped.
        dedup(bool): Remove repeated index columns (one per file).
    Returns:
        pandas.Dataframe
    """
    try:
        flist.aliases[0]
    except AttributeError:
        flist = FilesList(flist)
    # Determine how many columns each file has.
    numofcols = count_columns(flist, colSep=colSep)
    # Fetch and incrementally merge.
    result = None
    for f, (myfile, myalias) in flist.enum():
        # List the columns and remove the index one from among them.
        cols = [i for i in range(0,numofcols[f])]
        df = get_columns(FilesList(files=[myfile], aliases=[myalias]), cols=cols,
                        colSep=colSep, header=header, merge=False, index=index)[0]
        if (f == 0):
            result = df
        else:
            if saveHeader:
                # Extract headers, merge headers and tables separately, then put merged header on top of the merged table.
                hnew = pd.merge(left=result.iloc[[0]], right=df.iloc[[0]], sort=False,
                                left_index=True, right_index=True, how="outer")
                result = pd.concat([hnew, pd.merge(left=result.iloc[1:,:], right=df.iloc[1:,:], how=type, on=None,
                                                   left_index=True, right_index=True,
                                                   sort=False, suffixes=('','_'+ myalias))],
                                    axis=0, ignore_index=False, sort=False)
            else:
                result = pd.merge(left=result, right=df, how=type, on=None,
                                  left_index=True, right_index=True,
                                  sort=False, suffixes=('','_'+ myalias))
    # In addition to the new row_ID columns, the index column was kept for each table. Drop them as redundant.
    # If the index columns are not exact duplicates (due to gappy rows),
    # dedup_columns can be used afterwards on the merged file).
    if dedup:
        index_cols = [col for col in result.columns if '_|' + str(index) in col]
        result.drop(columns=index_cols, inplace=True)
    return result


# Helper function
def getDuplicateColumns(df):
    '''
    Get a list of duplicate columns.
    It will iterate over all the columns in dataframe and find the columns whose contents are duplicate.

    Stolen from https://thispointer.com/how-to-find-drop-duplicate-columns-in-a-dataframe-python-pandas/ .

    Args:
    	df: Dataframe object
    Returns:
        List of columns whose contents are redudnant (one occurence will of each will not be included in the list).
    '''
    duplicateColumnNames = set()
    # Iterate over all the columns in dataframe
    for x in range(df.shape[1]):
        # Select column at xth index.
        col = df.iloc[:, x]
        # Iterate over all the columns in DataFrame from (x+1)th index till end
        for y in range(x + 1, df.shape[1]):
            # Select column at yth index.
            otherCol = df.iloc[:, y]
            # Check if two columns at x 7 y index are equal
            if col.equals(otherCol):
                duplicateColumnNames.add(df.columns.values[y])
    return list(duplicateColumnNames)


def dedup_columns(flist, cols=[0,1], colSep=["\t"], merge=True):
    """Merge duplicate columns from the files, as they are.

    This function also supports key-aware appending, using outer-join, when a
    row index is specified.

    Args:
        flist: A list/FilesList of files to combine.
        colSep[str]: A list of characters used as field delimiters.
                    (Default ["\t"])
        cols[int] : A list of positional indexes of the
                    desired columns. (Default [0,1]).

    Returns:
        pandas.Dataframe
    """
    try:
        flist.aliases[0]
    except AttributeError:
        flist = FilesList(flist)
    # Determine how many columns each file has.
    numofcols = count_columns(flist, colSep=colSep)
    # Delegate fetching all the columns.
    result = []
    keyhead = None
    for f, (myfile, myalias) in flist.enum():
        # List the columns.
        allcols = [i for i in range(0,numofcols[f])]
        df = get_columns(FilesList(files=[myfile], aliases=[myalias]), cols=allcols,
        						colSep=colSep, header=False, merge=False, index=None)[0]
        # Collect duplicated values, drop duplicated columns, assign values to new column.
        v = df.iloc[:, cols].apply(lambda x: next(s for s in x.unique() if s), axis=1)
        df.drop(df.columns[cols], axis=1, inplace=True)
        df = pd.concat([df, v], axis=1, join='outer', sort=False, ignore_index=False)
        if not keyhead:
            keyhead = df.index.name
        result.append(df)
    # Merge.
    if merge:
        result = [pd.concat(result, axis=1, join='outer', ignore_index=False, sort=False), ]
        result[0].index.name = keyhead
    return result


def get_crosspoints(flist, cols=[0], rows=[0], colSep=["\t"], header=False, index=None, merge=True):
    """ Get the values at selected rows and columns.

    The values at the intersections of the selected rows and columns are extracted.

    Args:
        flist: A [str] list or fileutilities.FilesList of delimited text files.
        colSep[str]: List of column separators.
        cols[int]: List of columns.
        rows[int]: List of rows.
        header(bool): Whether there is a header line (False).
        index(int): Which column has the row labels (None).
        merge(bool): Merge results into single table (True).
    Returns:
        [pandas.DataFrame]:
    """
    results = get_columns(flist, cols=cols, colSep=colSep, header=header, merge=merge, index=index)
    for i in range(0, len(results)):
        results[i] = results[i].iloc[rows,:]
    return results


def store_metadata(flist, numoflines):
    """Store top lines of files into dictionary.

    Args:
        flist: A list or FilesList.
        numoflines(int): Number of lines to save.
    Returns:
        dict[]: The items of flist[] are used as keys.
    """
    metadata = dict()
    for myfile in flist:
        if myfile is None:
            fin = sys.stdin
        else:
            fin = open(myfile)
        lines = []
        for i in range(0, numoflines):
            lines.append(fin.readline())
        metadata[myfile] = "".join(lines)
        if fin != sys.stdin:
            fin.close()
    return metadata



#####   C L A S S E S   #####



class FilesList(list):
    """A container for a list of files.

    An extension of the built-in list, specifically for files, providing a
    means to import multiple filenames either from text lists or from
    directories. The purpose is to facilitate batch operations and sensible
    output of their results.

    FilesList is generally backwards compatible with the built-in list and it
    should be possible for them to be used interchangeably in most cases. A
    plain list can be cast as a FilesList, when necessary, allowing appointment
    of default alias values. A FilesList should always work as a plain list
    without additional actions (except for its string representation). When a
    FilesList is accessed as a plain list, only the full paths will be
    accessed. Certain equivalent methods are supplied for

    Most basic operations inherited from list are supported. Appending has been
    overridden to keep paths and aliases in sync. Sorting, deleting and
    inserting of items are currently not supported and will break the
    correspondence between paths and aliases.

    Attributes defined here:
        aliases = [] : Practical aliases for the full file-paths.
    """
    def __init__(self, files=None, aliases=None, fromtuples=None, verbatim=True):
        """Construct an instance of the FilesList.

        A FilesList can be created:
        - empty
        - from a list of files (with default aliases automatically assigned)
        - from a list of files and a list of aliases (in the same order)
        - from a list of (file, alias) tuples.

        Args:
            verbatim(bool): Whether to skip path pre-preocessing. (Default True)
            files[str]: A list of files. (Default None)
            aliases[str]: A list of aliases. (Default None)
            fromtuples[(str,str)]: A list of tuples (file, alias). (Default
                        None) If this is specified together with flist and/or
                        aliases, the data in fromtuples is used only.
        """
        # If data is a list of (file, alias) tuples, unpair tuples into lists.
        if fromtuples is not None:
            data = [list(t) for t in zip(*fromtuples)]
            # Any data passed to flist and aliases at method call is discarded.
            files = data[0]
            aliases = data[1]
        # Having aliases without the paths is rather useless.
        if aliases:
            if not files:
                raise ValueError("No files supplied for the aliases.")
        else:
            # Make None into empty.
            aliases = []
        # Assign default aliases to be same as files. Expand file paths.
        if files is not None:
            if not verbatim:
                files = expand_fpaths(files)
            if not aliases:
                for f in files:
                    aliases.append(self.autoalias(f))
        else:
            # If still empty, it was an empty call to the constructor.
            files = []
        # Create the basic list.
        super().__init__(files)
        # Add a plain list attribute for the aliases with default values.
        self.aliases = autonumerate(aliases)

    def __str__(self):
        """Represent as string.

        Overrides built-in list's representation.

        Returns:
            str
        """
        tmp = []
        for f, (myfile, myalias) in self.enum():
            tmp.append("\t".join([str(f), myfile, myalias]))
        tmp.append("")
        return "\n".join(tmp)

    def to_file(self, outfile=None, mode ='a'):
        """Save list as a text file that can be read back in.

        Args:
            outfile(str): Output file to write into. If omitted, it only
                        returns the content as a print-ready string.
                        (Default None)
            mode(str): Append ('a') or overwrite ('w'). (Default 'a')
        Returns:
            str: A print-ready multi-line string. This is returned even when an
                        output file is specified and written into.
        """
        result = ""
        for f, (myfile, myalias) in self.enum():
            result += myfile + "\t" + myalias + "\n"
        if outfile is not None:
            with open(outfile, mode) as out:
                out.write(result)
        return result

    def enum(self):
        """Enumerate as (index, (filepath, filealias)).

        Returns:
            enumerator"""
        return enumerate(zip(self, self.aliases))

    def get(self, loc):
        """Access path and alias at specified location as tuple.

        Args:
            loc[int]: Index of item to get.
        Returns:
            (str,str): Tuple (path, alias).
        """
        return (self[loc], self.aliases[loc])

    def append(self, myfile, myalias=None, verbatim=True):
        """Appends value to both the paths list and the aliases list.

        This method overrides the built-in append() of list. It is backwards
        compatible by automatically guessing an alias.
        This reduces the risk of the paths and aliases going out of sync due
        to items being manually added without updating the aliases.
        It is still possible to break the sync by manually adding items to the
        aliases.

        Args:
            myfile(str): File (path will be expanded).
            myalias(str): Alias for the file (Default None).
            verbatim(bool): Do not pre-process path for the target value. (Default True)
        """
        if myfile is not None:
            if not verbatim:
                myfile = expand_fpaths([myfile])[0]
        super().append(myfile)
        if not myalias:
            myalias = self.autoalias(myfile)
        self.aliases.append(myalias)
        self.aliases = autonumerate(self.aliases)

    def populate_from_files(self, myfiles, colSep="\t", verbatim=True):
        """Parse the list of files from one or multiple text files.

        Read in multiple lists of files from text and append them to the
        FilesList. All paths are automatically expanded and converted to
        absolute paths. Because of this, each file may also be given a more
        convenient alias. If no alias is given, the filename as supplied is
        used as the alias instead. The paths are not checked for existence.

        Existing contents of the object are kept and the new contents are
        appended.

        Input file format (no spaces allowed inside names):

            #comment
            path1/file1     alias1-prefix    alias1-suffix1     alias1-suffix2
            path1/file2     alias1-prefix    alias1-suffix3
            path2/file3     alias3
            path3/file4

        Args:
            file[str]: A list of text files each containing a list of files.
            colSep(str): Column separator. (Default "\\t")
            verbatim(bool): Do not pre-process paths for the target values. (Default True)
        Returns:
            FilesList: Returns self, to facilitate instantiation shortcuts.
        """
        # Read new list.
        paths = []
        for myfile in myfiles:
            with open(myfile, 'rU') as input:
                for line in input:
                    if line == "\n":
                        # Skip empty lines.
                        continue
                    elif line[0] == '#':
                        # Skip comments.
                        continue
                    else:
                        fields = line.rstrip().split(colSep)
                        paths.append(fields[0])
                        # Store the alias for the file.
                        if len(fields) > 1:
                            self.aliases.append("_".join(fields[1:]))
                        # If an alias was not specified, re-use the filepath given.
                        else:
                            self.aliases.append(self.autoalias(fields[0]))
        # Expand to absolute paths and add to main self list.
        if not verbatim:
            paths = expand_fpaths(paths)
        self.extend(paths)
        self.aliases = autonumerate(self.aliases)
        return self

    def populate_from_directories(self, dirpaths, patterns=None, verbatim=True):
        """Get files based on naming patterns from within a list of directories.

        Useful for selecting among files that follow a naming convention. The
        convention is represented by a list of regex strings, at least one of
        which has to match.
        File paths will be expanded automatically. The filenames will be used
        as the aliases.

        Existing contents of the object are kept and the new contents are
        appended.

        Args:
            dirpaths[str]: A list/FilesList of paths to directories from where
                        to get files.
            patterns[str]: A list of regex strings. Only files matching at least
                        one of these will be returned. The patterns will be
                        matched anywhere in the filenames.
            verbatim(bool): Do not pre-process paths for the target values. (Default True)
        Returns:
            FilesList: Returns self, to facilitate instantiation shortcuts.
        """
        rx = []
        if patterns:
            rx = [re.compile(p) for p in patterns]
        for d in dirpaths:
            try:
                os.chdir(d)
                for f in os.listdir(d):
                    if f in ["","\n",".",".."]:
                        continue
                    if not patterns:
                        # No filter.
                        self.append(os.path.join(d, f), self.autoalias(f), verbatim=verbatim)
                    else:
                        for p in rx:
                            if p.search(f):
                                self.append(os.path.join(d, f), self.autoalias(f), verbatim=verbatim)
                                break
            finally:
                pass
        self.aliases = autonumerate(self.aliases)
        return self.sorted()

    # Helper function.
    @staticmethod
    def autoalias(pathname):
        """Strip a path to the base filename."""
        if pathname is None:
            return None
        else:
            return os.path.splitext(os.path.basename(str(pathname)))[0]

    def sorted(self):
        """Sorted copy.

        Returns:
            FilesList
        """
        d = dict()
        for i, (myfile, myalias) in self.enum():
            d[myfile] = myalias
        sk = natural_sorted(list(d.keys()))
        newFL = FilesList()
        for k in sk:
            newFL.append(k, d[k])
        return newFL




#####   M A I N   #####


def main(args):
    """Provide command-line access to the module's functionality.

    The functionality and format of main is subject to change, as the module
    expands and evolves to suit my needs. Main() is not intended to be called
    from within other code.

    Optional short info is printed in commented lines. The info always
    *succeeds* the relevant output, rather than precede it. This serves as
    confirmation that the task completed. Calling details of the script are
    recorded in commented lines at the top of the output.

    For more info on the functionality, read the above documentation of the
    classes and functions. For usage syntax, execute the module with the -h
    argument.


    """
    # Organize arguments and usage help:
    parser = argparse.ArgumentParser(description="Provide INPUTTYPE and TARGETs \
     *before* providing any of the other parameters. This is due to many \
    parameters accepting an indefinite number of values. Only one task at a time.")

    # Input/Output.
    parser.add_argument('INPUTTYPE', type=str, choices=['L','T','D','P'],
                                help=" Specify the type of the TARGETs: \
                                'T' = The actual input filess. \
                                'L' = Text file(s) listing the input files. \
                                'P' = Get list of input files from STDIN pipe. \
                                'D' = Input data directly from STDIN pipe. \
                                ('D' is compatible with only some of the functions)")
    parser.add_argument('TARGET', type=str, nargs='*',
                                help=" The targets, space- or comma-separated. Usually files. \
                                Look into the specific task details below for special uses. \
                                Do not specify with INPUTTYPE 'P' or 'D'.")
    parser.add_argument('-O','--out', type=str, nargs=3,
                                help=" Send individual outputs to individual files instead of \
                                merging them to STDOUT. Output files will be like \
                                <out[0]>/<out[1]>target<out[2]>, where target is stripped of \
                                any directory path and its outermost file extension.")
    # Parameters.
    parser.add_argument('-L','--log', action='store_true',
                                help=" Log this command to ./commands.log.")
    parser.add_argument('-c','--comments', action='store_true',
                                help=" Include commented info to STDOUT or files. (Default don't include)")
    parser.add_argument('-C','--STDERRcomments', action="store_false",
                                help=" Do NOt show info in STDERR. (Default show)")
    parser.add_argument('-s','--sep', type=str, default=["\t"], nargs='+',
                                help=" A list of input field separators. The first value \
                                will be used for all output. (Default \\t, bash syntax for tab: $'\\t').")
    parser.add_argument('-l','--labels', action='store_true',
                                help=" Discard column headers (first content line) in input files. (Default do not discard)")
    parser.add_argument('-r','--relabel', action='store_false',
                                help=" Do not create new column headers that reflect the origin of the columns. (Default create)")
    parser.add_argument('-i','--index', action='store_true',
                                help=" Use column 0 as row index. The index will always be included in the output. (Default no index)")
    parser.add_argument('-M','--metadata', type=int, default=0,
                                help=" Number of metadata lines at the \
                                beginning of input data (Default 0). Metadate will be read separately \
                                and re-added verbatim into the output.")
    parser.add_argument('-R','--expand_ranges', action='store_true',
                                help=" If numeric ranges exist among the targets expand them as individual vlaues. \
                                Ranges must be in from:to format, inclusive of both end values. (Default False)")
    parser.add_argument('-V','--verbatim', action='store_true',
                                help=" Preserve the target values from a list file, do not try to expand them into absolute paths. (Default impute absolute paths)")
    # General tasks.
    parser.add_argument('--dir', type=str, nargs='*',
                                help=" List the contents of the target paths. \
                                Full absolute file paths are returned. Each file is also given an alias. \
                                Supplying an optional list of regex patterns enables filtering of the result.")
    parser.add_argument('--link', type=str, nargs='+',
                                help=" Create symbolic links for the targets into the specified directory. \
                                Any additional values are used as respective names for the links, one for one, \
                                otherwise the aliases or basenames will be used, enumerated when necessary.")
    parser.add_argument('--loop', type=str, nargs='+',
                                help=" Repeat the specified shell command for each target value. \
                                Available PLACEHOLDERS to insert the targets into the commands: \
                                {abs} full path, {dir} path of directory portion, {val} target value such as filename, \
                                {bas} basename (filename minus outermost extension), {ali} file alias. \
                                Flags intended for the nested command should be preceded \
                                by a ',' sign like this: ',-v'. Recursive calls to fileutilities.py are possible by \
                                nesting the placeholders and escapes: i.e. {{abs}}, ,,-v. One layer is peeled off \
                                with each call to fileutilities loop. The placeholders will take the values \
                                of the targets of the respectively nested call.")
    # Delimited file tasks.
    parser.add_argument('--concat', type=str,
                                help="Create an X-separated list out of the target values, where X is the string specified as argument here. \
                                Useful for creating comma-separated lists of files.")
    parser.add_argument('--swap', type=str,
                                help=" Replace all occurrences of the --sep values with the value supplied here.\
                                ** Bash syntax for tab: $'\\t'. Compatible with 'D' as INPUTTYPE.")
    parser.add_argument('--cntcols', action='store_true',
                                help="Count the number of fields in the first row of each target file.")
    parser.add_argument('--cols', nargs='+',
                                help="Extract the specified columns (named or 0-indexed) from each target. \
                                Column ranges in x:y format closed at both ends. \
                                Negative indices must be escaped first: \-1. Compatible with 'D' as INPUTTYPE.")
    parser.add_argument('--rndcols', type=int,
                                help="Randomly select this many columns from the target files. \
                                With --index, the index column will not be part of the random selection.")
    parser.add_argument('--appnd', type=str,
                                help="Append all the columns from same-length target files into a single table. \
                                Can be 'outer' or 'inner' join. If index is used, the values must be unique \
                                within each file.")
    parser.add_argument('--merge', nargs=3, type=str,
                                help="Merge table files. Index may contain duplicates and files may differ in dimensions. \
                                The first column of each file will be used as row index to merge on regardless of -i flag. \
                                First argument is type: 'left', 'right', 'inner', 'outer'. \
								Second argument is preserve first row: 'yes', 'no' (because merge sorts rows) \
								Third argument is deduplicate: 'yes', 'no' (index columns get repeated for each file. \
                                For very large tables, best don't deduplicate, and use --mrgdups on the output). \
                                For left/right joins, the order of files affects the result.")
    parser.add_argument('--mrgdups', type=int, nargs='+',
    							help="Combine gappy duplicate columns into a single column with all the values.\
    							Columns are specified by their 0-based positional index given as arguments here.")
    parser.add_argument('--valset', nargs=3,
                                help="Get the non-redundant set of values in the given row/column. \
                                Takes three arguments: (i) orientation 'r' for row or 'c' for column, \
                                (ii) position index of the row/column, (iii) repetition filter: \
                                'a' all values, 'u' unique values only, 'r' only values with two or more instances.")
    params = parser.parse_args(args)

    # INPUT ###################################################################

    targets = []
    for t in params.TARGET:
        v = t.split(",")
        if len(v) == 1:
            targets.append(t)
        else:
            targets.extend(v)

    flist = None
    if params.INPUTTYPE == 'P':
        # Read files list from STDIN
        flist = FilesList()
        for line in sys.stdin:
            fields = line.rstrip("\n").split("\t")
            if fields[0] != "":
                try:
                    flist.append(fields[0], fields[1])
                except IndexError:
                    flist.append(fields[0])
    elif params.INPUTTYPE == 'L':
        # Create the FilesList, by appending the contents of all provided lists.
        flist = FilesList().populate_from_files(targets, verbatim=params.verbatim)
    elif params.INPUTTYPE == 'T':
        # Create the FilesList by supplying a direct list of files.
        flist = FilesList(targets, verbatim=params.verbatim)
    elif params.INPUTTYPE == 'D':
        # Data will be read from STDIN. No files needed. Make an empty list.
        # Not all functions will switch to STDIN given this. Several will simply do nothing.
        flist = FilesList(verbatim=params.verbatim)
    else:
        sys.exit(ml.errstring("Unknown INPUTTYPE."))

    if params.expand_ranges:
        # Generate the range.
        myrange = []
        for t in params.TARGET:   # Look for multiple ranges.
            v = t.split(":")
            if len(v) > 1:
                myrange.extend(list(range(int(v[0]), int(v[1]) + 1)))
            else:
                myrange.extend(t)   # Assume string literal. Allows mixing numeric and string target values.
        flist = FilesList(myrange, verbatim=True)  # If some values are mere numbers, it's unlikely that any others are file paths.

    # Metadata. ---------------------------------------------------------------
    metadata = ""
    if params.metadata:
        metadata = store_metadata(flist, params.metadata)

    # OUTPUT ##################################################################

    outdir, outpref, outsuff = None, None, None
    if params.out:
        outdir = expand_fpaths([params.out[0]])[0]
        outpref = params.out[1]
        outsuff = params.out[2]

    # CALL DETAILS ############################################################

    if params.log:
        ml.log_command()
#    if params.STDERRcomments:
#        sys.stderr.write(ml.paramstring())

    # TASKS ###################################################################

    # Filter DIRECTORY contents. ----------------------------------------------
    if params.dir is not None:
        result = FilesList().populate_from_directories(flist, params.dir)
        try:
            if params.comments:
                sys.stdout.write(ml.paramstring())
            sys.stdout.write(result.to_file())
            if params.STDERRcomments:
                sys.stderr.write(ml.donestring("listing"))
        except IOError:
            pass


    # LOOP arbitrary command. -------------------------------------------------
    elif params.loop:
        command = []
        for c in params.loop:
            command.append(c.lstrip("+"))
        try:
            do_foreach(flist, command, out=(outdir, outpref, outsuff),
                       progress=(params.STDERRcomments), comments=params.comments,
                       log=params.log)
            if params.STDERRcomments:
                sys.stderr.write(ml.donestring("looping"))
        except IOError:
            pass


    # Symbolic LINKS. ---------------------------------------------------------
    elif params.link:
        slink(flist, dir=params.link[0], aliases=params.link[1:])
        if params.STDERRcomments:
            sys.stderr.write(ml.donestring("linking"))


    # CONCATENATE strings. --------------------------------------------------------
    if params.concat:
        sys.stdout.write(params.concat.join(flist))
        if params.STDERRcomments:
            sys.stderr.write(ml.donestring("concatenating values"))


    # SWAP substrings. --------------------------------------------------------
    elif params.swap is not None:
        result = swap_strFiles(flist, insep=params.sep, outsep=params.swap)
        # Create output filenames, if applicable. If [], then STDOUT.
        outfiles = make_names(flist.aliases, (outdir, outpref, outsuff))
        outstream = sys.stdout
        # I need the for loop to iterate at least once. Relevant for STDIN input, since I have no input files listed then.
        if flist == []:
            flist.append("<STDIN>")
        # Print the converted data.
        for i, (myfile, myalias) in flist.enum():
            if outfiles:
                # Send to individual file instead of STDOUT.
                outstream = open(outfiles[i], 'w')
            try:
                if params.comments:
                    # Embed call info at beginning of output. More useful there when outputting to files.
                    outstream.write(ml.paramstring("SOURCE: " + myfile))
                outstream.write(result[i].rstrip("\n") +"\n")
            except IOError:
                pass
            finally:
                if outfiles:
                    # Don't want to accidentally close STDOUT.
                    outstream.close()
        if params.STDERRcomments:
            try:
                sys.stderr.write(ml.donestring("swapping delimiters"))
            except IOError:
                pass


    # Get COLUMNS or RANDOM columns. (most code shared) -----------------------
    elif params.cols or params.rndcols:
        # Create output filenames, if applicable. If [], then STDOUT.
        outfiles = make_names(flist.aliases, (outdir, outpref, outsuff))
        outstream = sys.stdout
        merge = False if outfiles else True
        # Determine if using index, and assign appropriate value.
        idx = None
        if params.index:
            idx = 0
        else:
            idx = None
        # Extract data.
        result = None
        if params.cols:
            cols = []
            for p in params.cols:   # space separated arguments
                cols.extend(p.split(","))  # comma separated arguments
            # Get the specified columns.
            result = get_columns(flist, cols=cols, colSep=params.sep,
                                header=params.labels, merge=merge, index=idx)
        else:
            # Get random columns.
            result = get_random_columns(flist, k=params.rndcols, colSep=params.sep,
                                        header=params.labels, merge=merge, index=idx)
        # I need the for loop to iterate at least once. Relevant for STDIN input, since I have no input files listed then.
        if flist == []:
            flist.append("<STDIN>")
        if merge:
            try:
                if params.comments:
                    # Embed call info at beginning of output.
                    outstream.write(ml.paramstring("SOURCE: " + myfile))
                if params.metadata:
                    # Dump all the metadata from all the merged input sources.
                    for i, (myfile, myalias) in flist.enum():
                        outstream.write(metadata[myfile])
                outstream.write( result[0].to_csv(header=params.relabel, index=params.index, sep=params.sep[0]))
            except IOError:
                pass
        else:
            for i, (myfile, myalias) in flist.enum():
                outstream = open(outfiles[i], 'w')
                try:
                    if params.comments:
                        # Embed call info at beginning of output.
                        outstream.write(ml.paramstring("SOURCE: " + myfile))
                    if params.metadata:
                        outstream.write(metadata[myfile])
                    outstream.write( result[i].to_csv(header=params.relabel, index=params.index, sep=params.sep[0]))
                except IOError:
                    pass
                finally:
                    outstream.close()
        if params.STDERRcomments:
            try:
                if params.cols:
                    sys.stderr.write(ml.donestring("getting columns, index "+ str(idx is not None)))
                else:
                    sys.stderr.write(ml.donestring("getting random columns, index "+ str(idx is not None)))
            except IOError:
                pass


    # APPEND columns or MERGE table. ---------------------------------------------------------
    elif params.appnd or params.merge:
        idx = None
        if params.appnd:
            if params.index:
                idx = 0
            df = append_columns(flist, colSep=params.sep, header=params.labels, index=idx, type=params.appnd)
        else:
            df = merge_tables(flist, colSep=params.sep, header=params.labels, index=0, type=params.merge[0],
                              saveHeader=params.merge[1] == "yes", dedup=params.merge[2] == "yes")
        try:
            if params.comments:
                ml.parastring()
            if params.metadata:
                # Dump all the metadata from all the merged input sources.
                for i, (myfile, myalias) in flist.enum():
                    outstream.write(metadata[myfile])
            sys.stdout.write(df.to_csv(sep=params.sep[0], header=params.relabel, index=params.index))
            if params.STDERRcomments:
                if params.appnd:
                    sys.stderr.write(ml.donestring(params.appnd +" append of columns, index "+ str(idx is not None)))
                else:
                    sys.stderr.write(ml.donestring(params.merge[0] +" merge of tables"))
        except IOError:
            pass


    # MERGE duplicate columns. ---------------------------------------------------------
    elif params.mrgdups:
        # Create output filenames, if applicable. If [], then STDOUT.
        outfiles = make_names(flist.aliases, (outdir, outpref, outsuff))
        outstream = sys.stdout
        merge = False if outfiles else True
        # Do.
        result = dedup_columns(flist, colSep=params.sep, cols=params.mrgdups)
		# I need the for loop to iterate at least once. Relevant for STDIN input, since I have no input files listed then.
        if flist == []:
            flist.append("<STDIN>")
        if merge:
            try:
                if params.comments:
                    # Embed call info at beginning of output.
                    outstream.write(ml.paramstring("SOURCE: " + myfile))
                if params.metadata:
                    # Dump all the metadata from all the merged input sources.
                    for i, (myfile, myalias) in flist.enum():
                        outstream.write(metadata[myfile])
                outstream.write( result[0].to_csv(header=params.relabel, index=params.index, sep=params.sep[0]))
            except IOError:
                pass
        else:
            for i, (myfile, myalias) in flist.enum():
                outstream = open(outfiles[i], 'w')
                try:
                    if params.comments:
                        # Embed call info at beginning of output.
                        outstream.write(ml.paramstring("SOURCE: " + myfile))
                    if params.metadata:
                        outstream.write(metadata[myfile])
                    outstream.write( result[i].to_csv(header=params.relabel, index=params.index, sep=params.sep[0]))
                except IOError:
                    pass
                finally:
                    outstream.close()
        if params.STDERRcomments:
            try:
                sys.stderr.write(ml.donestring("deduplicating columns"))
            except IOError:
                pass


    # COUNT columns. ----------------------------------------------------------
    elif params.cntcols:
        result = count_columns(flist, params.sep)
        try:
            if params.comments:
                sys.stdout.write(ml.paramstring())
            for f, (myfile, myalias) in flist.enum():
                print("\t".join([str(result[f]), myalias, myfile]))
            if params.STDERRcomments:
                sys.stderr.write(ml.donestring("counting columns"))
        except IOError:
            pass


    # SET of values in row/column. --------------------------------------------
    elif params.valset:
        nest = get_valuesSet(flist, axis=params.valset[0], index=params.valset[1], filter=params.valset[2], colSep=params.sep)
        try:
            if params.comments:
                sys.stdout.write(ml.paramstring())
            for f, (myfile, myalias) in flist.enum():
                print("".join([myfile, "\t", str(nest[f])]))
            if params.STDERRcomments:
                sys.stderr.write(ml.donestring("obtaining set of values."))
        except IOError:
            pass



#     # All done.
#     if params.STDERRcomments:
#         sys.stderr.write(ml.donestring())





#####    E X E C U T I O N   #####


# Call main only if the module was executed directly.
if __name__ == "__main__":
    main(sys.argv[1:])

    sys.exit(0)

#EOF
